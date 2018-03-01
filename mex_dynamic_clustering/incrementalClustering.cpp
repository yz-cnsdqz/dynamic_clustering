#include <algorithm>
#include <chrono>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <sstream>
#include <string>
#include <mutex>
#include <string>
#include <vector>
#include <cmath>
#include <utility> //std::pair

#include <pthread.h>
#include <time.h>
#include <signal.h>
#include <stdio.h>  // snprintf
#include <unistd.h>
#include <stdlib.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <sys/types.h>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/array.hpp>
#include <gflags/gflags.h>

#include <opencv2/core/core.hpp>
// #include <opencv2/contrib/contrib.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "Persistence1D/src/persistence1d/persistence1d.hpp"
#include "mex.h"
// #include "matrix.h"
#include "opencvmex.hpp"


// Global parameters
// int DISPLAY_RESOLUTION_WIDTH;
// int DISPLAY_RESOLUTION_HEIGHT;
// int CAMERA_FRAME_WIDTH;
// int CAMERA_FRAME_HEIGHT;
// int NET_RESOLUTION_WIDTH;
// int NET_RESOLUTION_HEIGHT;
// int BATCH_SIZE;
// float SCALE_GAP;
// float START_SCALE;
// int NUM_GPU;
// std::string PERSON_DETECTOR_CAFFEMODEL; //person detector
// std::string PERSON_DETECTOR_PROTO;      //person detector
// std::string POSE_ESTIMATOR_PROTO;       //pose estimator
// const auto MAX_PEOPLE = RENDER_MAX_PEOPLE;  // defined in render_functions.hpp
// const auto BOX_SIZE = 368;
// const auto BUFFER_SIZE = 4;    //affects latency
// const auto MAX_NUM_PARTS = 70;


      
// global queues for I/O
struct Global {

    int num_samples = 0;
    int num_dims = 0;
    const float thresh_from_buffer = 0.5;
    int KEYFRAME_BUFFERSIZE = 0;       // BUFFER size. Reaching this size triggers [ALGO: keyframe selection] and then release the buffer.
    static const int KEYFRAME_SAMPLING_RATE = 1;
    int n_frames_inBuffer = 0;  // after processing the frames in the buffer, the first X frames are erased from the vector
    // int n_frames_inPool = 0;
    static const int MIN_KEYFRAME_POOLSIZE = 10;   //  Achieving this number triggers clustering algorithm and then releases this buffer. 
    int KEYFRAME_POOL_FILLTIMES = 0;    // when the keyframe pool is full, then perform clustering and then released, this number +1.
    static const bool USE_STANDARDIDATION = false;
    int USE_TEMPORAL_REG = 0;
    float temporal_reg_alpha = 0;
    int VERBOSE = 0;
    static const bool SAVE_IMAGE = true;
    static const bool SAVE_CODEBOOK=true;
    std::vector<cv::Mat> output_vec_KFrameBuffer; // a set of frames are stored within the buffer ( not always from the beginning)
    std::vector<cv::Mat> output_vec_KFramePool;   // the set of keyframe candidiates. The sorted result
    static const bool IS_DEBUG = false;
    float RADIUS_OFFSET = 0; // compensate floating error
    struct Codebook {
        int n_clusters;
        cv::Mat init_standardization_mean;
        cv::Mat init_standardization_std;
        cv::Mat init_sample_labels; 
        cv::Mat sample_labels;
        const int transition_map_dims = 2;
        const float init_cluster_std = 0.0f;
        int init_k_for_kmeans = 0;
        std::vector<int> n_samples_per_cluster;
        
        cv::Mat cluster_labels;
        // std::vector<Frame> cluster_imgs = std::vector<Frame> (); // the closest samples to the centers
        cv::Mat cluster_locs; // we can append row vectors later
        cv::Mat cluster_std; // This vector has the same length of n_clusters. Each element equals to 3*max(sigma)
        
        cv::SparseMat transition_map;
        int index_cluster_pre;
        cv::Mat prob_prediction;

    };
    Codebook codebook;
 };



template<typename T>
void test_fun_Cout(T& a){
    std::cout << a << std::endl;
}

struct ColumnCompare
{
    bool operator()(const std::vector<float>& lhs,
                    const std::vector<float>& rhs) const
    {
        return lhs[2] > rhs[2];
        //return lhs[0] > rhs[0];
    }
};

Global global;



float get_wall_time() {
    struct timeval time;
    if (gettimeofday(&time,NULL)) {
        //  Handle error
        return 0;
    }
    return (float)time.tv_sec + (float)time.tv_usec * 1e-6;
    //return (float)time.tv_usec;
}



template<typename T>
class SortingIndexExtractor{
std::vector<T> _vec;
public:
    SortingIndexExtractor(std::vector<T>& vec_in) : _vec(vec_in){}
    bool operator () (const int& a, const int& b) const {
        return _vec[a] < _vec[b];
    }
};



void calClusterRadiusfromStd( std::vector<float>& radius ) {
    for (int i = 0; i < global.codebook.n_clusters; i++) {
        cv::Mat _std = global.codebook.cluster_std.row(i);
        double _std_min, _std_max;
        cv::minMaxIdx(_std, &_std_min, &_std_max);
        radius.push_back((float)_std_max * 3.0f + global.RADIUS_OFFSET);
    }
}

void calClusterStd(cv::Mat& samples) {
    // std::cout<< "[DEBUG] calClusterRadius(cv::Mat& samples) starts..." <<std::endl;
    global.codebook.n_samples_per_cluster = std::vector<int> (global.codebook.n_clusters, 0); 
    global.codebook.cluster_std = cv::Mat(global.codebook.cluster_locs.size(), CV_32F);
    // std::cout<< "[DEBUG] calClusterRadius(cv::Mat& samples) loop starts..." <<std::endl;
    for (int i = 0; i < global.codebook.n_clusters; i++){
        cv::Mat _samplei;       

        for (int j = 0; j < samples.rows; j++) {
            // std::cout << "[DEBUG] global.codebook.init_sample_labels.at<int>(j) = " << global.codebook.init_sample_labels.at<float>(j) <<std::endl;
            if (global.codebook.init_sample_labels.at<float>(j)==(float)i){
                _samplei.push_back(samples.row(j));    
                global.codebook.n_samples_per_cluster[i] +=1;
            }

        }
        // std::cout<< "[DEBUG] calClusterRadius(cv::Mat& samples) global.codebook.n_samples_per_cluster[0] = " 
        //          << global.codebook.n_samples_per_cluster[0] <<std::endl;
        int n_samplei = global.codebook.n_samples_per_cluster[i];
        // std::cout << n_samplei << std::endl;
        // std::cout<< "[DEBUG] calClusterStd(cv::Mat& samples) _samplei.size() = " << _samplei.size() <<std::endl;
        cv::Mat _mean = global.codebook.cluster_locs.row(i);
        cv::Mat _mean_repmat = cv::repeat(_mean, n_samplei, 1);
        // std::cout << _mean_repmat <<std::endl;
        cv::Mat _std;
        cv::reduce ((_samplei-_mean_repmat).mul((_samplei-_mean_repmat))/n_samplei, _std, 0, CV_REDUCE_SUM);
        // std::cout<< "[DEBUG] calClusterRadius(cv::Mat& samples) after reduce..." <<std::endl;
        cv::sqrt(_std,_std);
        _std = _std + global.codebook.init_cluster_std;
        _std.copyTo(global.codebook.cluster_std.row(i));
        
    }
    // std::cout<< "[DEBUG] calClusterRadius(cv::Mat& samples) ends..." <<std::endl;
}


void sampleStandardization(cv::Mat& samples, cv::Mat& samples_dst, bool flag) {
    
    if (flag){
    // 1. Calculate mean and std along each dimension
        cv::Mat _mean, _mean_repmat;
        cv::Mat _std, _std_repmat;
        cv::reduce(samples, _mean, 0, CV_REDUCE_AVG);
        // std::cout<< "[DEBUG] _mean = " << _mean <<std::endl;
        cv::repeat(_mean, samples.rows, 1, _mean_repmat);
        cv::reduce((samples-_mean_repmat).mul(samples-_mean_repmat), _std, 0, CV_REDUCE_AVG);
        cv::sqrt(_std, _std);
        cv::repeat(_std, samples.rows, 1, _std_repmat);
        global.codebook.init_standardization_mean = _mean.clone();
        global.codebook.init_standardization_std = _std.clone();
        
        // 2. perform standardization
        cv::divide(samples-_mean_repmat, _std_repmat, samples_dst);
    }
    else{
        samples_dst = samples.clone();
        global.codebook.init_standardization_mean = cv::Mat::zeros(1,global.num_dims, CV_32F);
        global.codebook.init_standardization_std = cv::Mat::ones(1,global.num_dims,CV_32F);
    }
    
}

// void* calClusterImgs( ) {
//     // std::cout << "imshowing starting" << std::endl;
//     // std::cout << "imshowing starting 2" << std::endl;
//     // std::cout << "[DEBUG] calClusterImgs() begins" <<std::endl;
//     cv::Mat samples_one_cluster;
//     // std::cout << "[debug] global.output_vec_KFramePool.size() = " << global.output_vec_KFramePool.size() << std::endl;

//     float _distn; 
//     int idx = 0;
//     int n_samples = global.KEYFRAME_POOLSIZE;
//     // std::cout << "n_samples = " << n_samples << std::endl;
//     // std::cout << "--break1--" <<std::endl;
    
//     // std::cout << "--break2--" <<std::endl;
//     for (int i = 0; i < global.codebook.n_clusters; i++){
//         // std::cout << "--cluster = " << i << std::endl;
//         cv::Mat c_loc = global.codebook.cluster_locs.row(i);
//         float _dist = 1e10;
//         for (int j = 0; j < n_samples; j++) {
//             // std::cout << "j=" <<j<<std::endl;
//             // std::cout << "label="<<global.codebook.cluster_labels.at<int>(j)<<std::endl;
//             if (global.codebook.init_sample_labels.at<int>(j)==i) {
//                 cv::Mat jts, jtst;
//                 // jts = global.output_vec_KFramePoolSkeleton[j];
//                 cvtFrameToJoints(global.output_vec_KFramePool[j], jts);
//                 jtst = jts.t();
//                 // std::cout << "count" <<std::endl;

//                 _distn = cv::norm(jtst, c_loc);
//                 // std::cout << _distn << std::endl;
//                 if ( _distn < _dist ) {
//                     _dist = _distn;
//                     idx = j;
//                 }
//             }
//         }

//         global.codebook.cluster_imgs.push_back(global.output_vec_KFramePool[idx]);
//         // std::cout << "[DEBUG] current center index=" << i << std::endl;   
//         // std::cout << "[DEBUG] original index in the pool=" << idx << std::endl;
//     }

//     // cv::namedWindow("codebook image");
//     // for (int k = 0; k < global.codebook.cluster_imgs.size(); k++) {
//     //     std::cout << "cluster - " << k+1 << std::endl;
//     //     cv::imshow("codebook image", global.codebook.cluster_imgs[k].img);
//     //     cv::waitKey(0);
//     // }

//     return nullptr;
// }

void showCodebookInfo( ) {
    std::cout << "------------------------Final Result ---------------------------" <<std::endl;
    std::cout << "[CODEBOOK INFO] codebook.n_clusters = " << global.codebook.n_clusters << std::endl;
    // std::cout << "[CODEBOOK INFO] codebook.cluster_imgs.size() = " << global.codebook.cluster_imgs.size() << std::endl;
    // FrameCompare comp;
    // std::priority_queue<Frame, std::vector<Frame>, FrameCompare> buffer(comp);
    // for (int j = 0; j < global.codebook.cluster_imgs.size(); j++) 
    //     buffer.push(global.codebook.cluster_imgs[j]);

    
    std::stringstream out_name;

    // cv::namedWindow("Codebook Image");
    for (int k = 0; k < global.codebook.cluster_locs.rows; k++) {
        std::cout << "[CODEBOOK INFO] cluster - " << k << std::endl;
        std::cout << "[CODEBOOK INFO] codebook.cluster_locs = " << global.codebook.cluster_locs.row(k) << std::endl;
        std::cout << "[CODEBOOK INFO] codebook.cluster_std = " << global.codebook.cluster_std.row(k)<< std::endl;    // std::cout << "imshowing starting 2" << std::endl;
        std::cout << "[CODEBOOK INFO] codebook.num_samples = " << global.codebook.n_samples_per_cluster[k]<< std::endl;    // std::cout << "imshowing starting 2" << std::endl;
        // cv::imshow("Codebook Image", global.codebook.cluster_imgs[k].img);
        // std::cout << "video_frame_index = " << global.codebook.cluster_imgs[k].video_frame_number << std::endl;
        // if(global.SAVE_IMAGE) {
        //     out_name << "output/codebook_cluster_" << k << ".png";
        //     std::cout << out_name.str()<<std::endl;
        //     cv::imwrite(out_name.str(), global.codebook.cluster_imgs[k].img);
        //     out_name.str("");
        // }
       
    }
}


/*
void* showBufferContents( ) {

    // while (1) {
    //     if (global.quit_threads)
    //         break;
    
    // lock the global memory, so other threads should wait for the following computation.    
    // std::unique_lock<std::mutex> lock{global.mutex};
    int n_frames = global.KEYFRAME_POOLSIZE;
    std::cout << "global.output_vec_KFramePool.size()=" << n_frames <<std::endl;
    std::stringstream out_name;    
    // sort the frames according to indeces in the ascending order.
    FrameCompare comp;
    std::priority_queue<Frame, std::vector<Frame>, FrameCompare> buffer(comp);
    for (int j = 0; j < n_frames; j++) 
        buffer.push(global.output_vec_KFramePool[j]);

    // show the results from the priority_queue
    cv::namedWindow("buffer content", CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
    for (int i = 0; i < n_frames; i++) {
        // std::cout << "global.output_vec_KFrameBuffer[end].index=" << global.output_vec_KFrameBuffer[i].index <<std::endl;
        // std::cout << "global.output_vec_KFrameBuffer[end].img.rows=" << global.output_vec_KFrameBuffer[i].img.rows <<std::endl;
        // std::cout << "global.output_vec_KFrameBuffer[end].img.cols=" << global.output_vec_KFrameBuffer[i].img.cols <<std::endl;
        // std::cout << "global.output_vec_KFrameBuffer[end].img.size().height=" << global.output_vec_KFrameBuffer[i].img.size().height <<std::endl;
        // std::cout << "global.output_vec_KFrameBuffer[end].img.size().width=" << global.output_vec_KFrameBuffer[i].img.size().width <<std::endl;
        Frame frame = buffer.top();


        redrawSkeletonsUsingJoint(frame);
        // cv::imshow("buffer content", frame.img);
        // test_fun_showMatData(frame);


        int offset = DISPLAY_RESOLUTION_WIDTH * DISPLAY_RESOLUTION_HEIGHT;
        for(int c = 0; c < 3; c++) {
            for(int i = 0; i < DISPLAY_RESOLUTION_HEIGHT; i++) {
                for(int j = 0; j < DISPLAY_RESOLUTION_WIDTH; j++) {
                    int value = int(frame.data_for_mat_update[c*offset + i*DISPLAY_RESOLUTION_WIDTH + j] + 0.5);
                    value = value<0 ? 0 : (value > 255 ? 255 : value);
                    frame.data_for_wrap[3*(i*DISPLAY_RESOLUTION_WIDTH + j) + c] = (unsigned char)(value);
                }
            }
        }

        cv::Mat wrap_frame(DISPLAY_RESOLUTION_HEIGHT, DISPLAY_RESOLUTION_WIDTH, CV_8UC3, frame.data_for_wrap);

        // std::cout << " - display the " << i+1 << "-th frame." <<std::endl; 

        // cv::imshow("buffer content", wrap_frame);
        out_name << "output/kframe_" << i << ".png";
        std::cout << out_name.str()<<std::endl;
        cv::imwrite(out_name.str(), wrap_frame);
        out_name.str("");

        // cv::waitKey(0);
        buffer.pop();
    }
        
    // lock.unlock();
    return nullptr;
}*/




int calClusterNumbers(cv::Mat& samples){
    int n_cluster = 0;
    cv::Mat _W (cv::Size(samples.rows, samples.rows), CV_32F, cv::Scalar(0));
    cv::Mat _D;
    cv::Mat _L;
    cv::Mat _L_eigenval;
    cv::Mat _aux_mat;
    int n_clusters;
    
    for (int i = 0; i < _W.rows; i++) {
        for (int j = i+1; j < _W.cols; j++){
            cv::Mat _jts1 = samples.row(i);
            cv::Mat _jts2 = samples.row(j);
            _W.at<float> (i,j) = cv::norm(_jts1, _jts2, cv::NORM_L2);
        }
    } 

    _W = (_W + _W.t())/2.0 ;
    // std::cout << "[DEBUG] _W = " << _W <<std::endl;
    cv::reduce(_W, _D, 0, CV_REDUCE_SUM, CV_32F);

    // _L = cv::Mat::diag(_D)-_W; // unnormalized graph laplacian
    _L =  cv::Mat::eye(_W.size(), CV_32F) - cv::Mat::diag(_D).inv() * _W; // normalized graph laplacian. See the literature

    
    cv::eigen(_L, _L_eigenval);
    cv::sort(_L_eigenval, _L_eigenval, CV_SORT_EVERY_COLUMN+CV_SORT_ASCENDING);
    cv::sortIdx( _L_eigenval(cv::Range(1,_W.rows),cv::Range::all())-_L_eigenval(cv::Range(0,_W.rows-1),cv::Range::all()), 
        _aux_mat, CV_SORT_EVERY_COLUMN+CV_SORT_DESCENDING);
    n_clusters = _aux_mat.at<int>(0)+1;    

    return n_clusters;
}







void temporalSpline(cv::Mat& src, cv::Mat& dst){
    //Input: src, Nxglobal.num_dims matrix
    //Output: des, Nxglobal.num_dims matrix
    //Algorithm: spline along each column
    
    int n_dims = src.cols;
    int n_samples = src.rows;
    static const int bt = 1;
    static const float w = 1.0f;
    cv::Mat src_mirror;

    cv::copyMakeBorder(src,src_mirror, bt,bt,0,0,cv::BORDER_REFLECT_101);
    cv::Mat mask = cv::Mat::zeros(src_mirror.size(), CV_32F);
    
    static const float eps = 1e-6f;
    static const float alpha = 0.1f;
    
    cv::threshold(src_mirror,mask, 0.01, 1.0, 0 );
    int iterMax = 80;
    
    for (int i = 0; i < n_dims; i++){
        cv::Mat f = src_mirror.col(i);
        cv::Mat u = f.clone();
        cv::Mat a = mask.col(i);
        for(int it = 0; it < iterMax; it++){
            for(int t = bt; t < bt+n_samples; t++){
                float psi_prime = 1.0f/sqrtf( std::pow(u.at<float>(t)-f.at<float>(t), 2) + eps);
                float phi_prime = 1.0f/( std::pow(0.5*u.at<float>(t+1)-0.5*u.at<float>(t-1), 2) + eps);
                // u.at<float>(t) = ( (1.0f-a.at<float>(t))*f.at<float>(t) 
                // + (eps+a.at<float>(t))*(4.0f*u.at<float>(t-1)+4.0f*u.at<float>(t+1)-u.at<float>(t-2)-u.at<float>(t+2)))/(1+5.0f*a.at<float>(t)+6.0f*eps);
                u.at<float>(t) = (1.0f-w)*u.at<float>(t) 
                                 + w*(a.at<float>(t)*psi_prime*f.at<float>(t) + alpha*phi_prime*(u.at<float>(t-1)+u.at<float>(t+1)))
                                 /(a.at<float>(t)*psi_prime+2.0f*alpha*phi_prime);
            }

            u.at<float>(bt-1) = u.at<float>(bt);    
            u.at<float>(bt+n_samples) = u.at<float>(bt+n_samples-1);
            
       }
        // dst.col(i) = u.clone();
        cv::Mat u_o = u.rowRange(bt, bt+n_samples);
        u_o.copyTo(dst.col(i));
    }
}




void temporalSpline1D(cv::Mat& src, cv::Mat& dst){
    //Input: src, Nxglobal.num_dims matrix
    //Output: des, Nxglobal.num_dims matrix
    //Algorithm: spline along each column
    
    int n_samples = src.rows;
    static const int bt = 1;
    static const float w = 1.96f;
    cv::Mat f, u;

    cv::copyMakeBorder(src,f, bt,bt,0,0,cv::BORDER_REFLECT_101);
    // a = cv::Mat::zeros(f.size(), CV_32F);
    u = f.clone();
    // u = cv::Mat::zeros(f.size(), CV_32F);
    
    static const float lambda = 1e-5f;
    float alpha = global.temporal_reg_alpha;
    
    // cv::threshold(f,a, 0.01, 1.0, 0 );
    int iterMax = 1000;    
    
    for(int it = 0; it < iterMax; it++){
        for(int t = bt; t < bt+n_samples; t++){
            float psi_prime = 1.0f/sqrtf( std::pow(u.at<float>(t)-f.at<float>(t), 2) + 1e-8);
            // float phi_prime = 1.0f/( std::pow(0.5*u.at<float>(t+1)-0.5*u.at<float>(t-1), 2) / lambda + 1);
            float phi_prime = 1.0f/sqrtf( std::pow(0.5*u.at<float>(t+1)-0.5*u.at<float>(t-1), 2) + 1e-8);
            // u.at<float>(t) = ( (1.0f-a.at<float>(t))*f.at<float>(t) 
            // + (eps+a.at<float>(t))*(4.0f*u.at<float>(t-1)+4.0f*u.at<float>(t+1)-u.at<float>(t-2)-u.at<float>(t+2)))/(1+5.0f*a.at<float>(t)+6.0f*eps);
            u.at<float>(t) = (1.0f-w)*u.at<float>(t) 
                             + w*(psi_prime*f.at<float>(t) + alpha*phi_prime*(u.at<float>(t-1)+u.at<float>(t+1)))
                             /(psi_prime+2.0f*alpha*phi_prime);
        }

        u.at<float>(bt-1) = u.at<float>(bt);    
        u.at<float>(bt+n_samples) = u.at<float>(bt+n_samples-1);
    }
    // dst.col(i) = u.clone();
    cv::Mat u_o = u.rowRange(bt, bt+n_samples);
    u_o.copyTo(dst);
}








float getBoundingBoxAreaFromJoint(const cv::Mat& joint){
    float xmin = 10e10;
    float ymin = 10e10;
    float xmax = 10e-10;
    float ymax = 10e-10;
    float bb_area;
    int N = joint.cols;
    cv::Mat joint_ascending, joint_descending;
    cv::Mat joint_ascending_idx, joint_descending_idx;
    
    cv::sort(joint, joint_ascending, CV_SORT_EVERY_ROW+CV_SORT_ASCENDING);
    cv::sort(joint, joint_descending, CV_SORT_EVERY_ROW+CV_SORT_DESCENDING);
    cv::sort(joint, joint_ascending_idx, CV_SORT_EVERY_ROW+CV_SORT_ASCENDING);
    cv::sort(joint, joint_descending_idx, CV_SORT_EVERY_ROW+CV_SORT_DESCENDING);
    for(int i = 0; i < N; i++){
        if ( (joint_ascending_idx.at<int>(i)%2==0)&&( xmin > joint_ascending.at<float>(i)) )
            xmin = joint_ascending.at<float>(i);

        if ( (joint_ascending_idx.at<int>(i)%2==1)&&( ymin > joint_ascending.at<float>(i)) )
            ymin = joint_ascending.at<float>(i);

        if ( (joint_descending_idx.at<int>(i)%2==0)&&( xmax < joint_descending.at<float>(i)) )
            xmax = joint_descending.at<float>(i);

        if ( (joint_descending_idx.at<int>(i)%2==1)&&( ymax < joint_descending.at<float>(i)) )
            ymax = joint_descending.at<float>(i);
    }

    bb_area = (xmax-xmin)*(ymax-ymin);

    return bb_area;
       
}


void calPosePattern (const cv::Mat& joints, cv::Mat& pose_pattern, 
    cv::Mat& extension, cv::Mat& movement){


    int n_frames  = joints.rows;
    int n_dims = joints.cols; // n_dims should be global.num_dims
    pose_pattern = joints.clone();
    extension = cv::Mat::zeros(joints.rows, 1, CV_32F);
    movement = cv::Mat::zeros(joints.rows, 1, CV_32F);
    cv::Mat movement_mean = cv::Mat::zeros(joints.rows, 1, CV_32F);
    
    cv::Mat movement_score = cv::Mat::ones(1, joints.cols, CV_32F); // Large part movement has more influence on the movement score.
    
    
    // cv::Mat _move = cv::Mat::zeros(joints.size(), CV_32F);
    for(int t = 0; t < n_frames; t++){
        
        // calculate the body extension, measured by bounding box area
        extension.at<float>(t) = getBoundingBoxAreaFromJoint(joints.row(t));
        // calculate the motion energy
        if (t==0 || t==n_frames-1)
            continue;
        else
            for(int j = 0; j < n_dims; j++)
                movement_mean.at<float>(t) += (std::fabs(0.5f*joints.at<float>(t+1,j)-0.5f*joints.at<float>(t-1,j)))/n_dims;


        
    } 
    movement = movement_mean.clone();    

 }


void calBipolarFeature(cv::Mat& src, cv::Mat& dst){
    cv::Mat mean, mean_repmat;
    cv::reduce(src, mean, 0, CV_REDUCE_AVG);
    cv::repeat(mean, src.rows, 1, mean_repmat);

    cv::Mat src_cen = src-mean_repmat;
    cv::threshold(src_cen, dst, 0, 2, cv::THRESH_BINARY);
    dst = dst-1.0f;
}

void calBidirectionalAssociativeMemoryEnergy(cv::Mat& X, cv::Mat& extension, cv::Mat& movement, cv::Mat& mem_energy){
    
    cv::Mat W = X.t()*X-(X.rows)*cv::Mat::eye(global.num_dims,global.num_dims,CV_32F);
    cv::Mat WW = -(0.5*X*W*X.t()).diag(0);


    mem_energy = WW.clone();

}




void getKeyFrameIdxFromBuffer(const cv::Mat& memory, const cv::Mat& movement, std::vector<int>& frame_idx){
    // Input: memory, movement
    // Output: frame_idx - the frame index set, the temporal order can be ignored since we have a sorting scheme based on timestamp later.

    
    std::vector<float> _curve;
    for(int i = 0; i < memory.rows; i++){
        // _memory.push_back(memory.at<float>(i));
        // _motion.push_back(movement.at<float>(i));
        // _curve.push_back(_memory.back() * _motion.back()); 
        _curve.push_back(float(memory.at<float>(i) + movement.at<float>(i)));
        // frame_idx.push_back(i);
    }
    cv::Mat curve_mat (_curve);
    curve_mat.convertTo(curve_mat, CV_32F);
    cv::normalize(curve_mat, curve_mat, 0,1, cv::NORM_MINMAX);
    
    for(int i = 0; i < global.KEYFRAME_BUFFERSIZE; i+= global.KEYFRAME_SAMPLING_RATE){
        if ( curve_mat.at<float>(i) < global.thresh_from_buffer)
            frame_idx.push_back(i);
    }

    SortingIndexExtractor<float> idx_sorter(_curve);
    sort(frame_idx.begin(), frame_idx.end(),idx_sorter);

    // // determine K for Kmeans. This module only works at the first working memory for initialisation
    // // However, the value of K is not equivalent to the temporal segmentation. Therefore, we create a fully connec
    if (global.KEYFRAME_POOL_FILLTIMES==0){
        
        // extract maxima in the memory+movement curve, so as to segment the video.
        p1d::Persistence1D p;
        p.RunPersistence(_curve);
        std::vector<p1d::TPairedExtrema> extrema;
        extrema.clear();
        p.GetPairedExtrema(extrema, 0.05);
        std::vector<int> segment_idx = {0};
        segment_idx.clear();
        for(std::vector<p1d::TPairedExtrema>::iterator it = extrema.begin(); it != extrema.end(); it++){
            segment_idx.push_back((*it).MaxIndex);
            std::cout << segment_idx.back() << std::endl;
        }
        global.codebook.init_k_for_kmeans = 0;
        int KK = segment_idx.size();
        // If no frames are extracted from a segment, then give up this segment.
        // Note that here we only consider the number of K. For evaluating the segmentation, one can check another git branch.
        int meter = 0;
        if (KK!=0){
            for(int i = 0; i < segment_idx.size()-1; i++){
                int a1 = segment_idx[i];
                int a2 = segment_idx[i+1];
                
                for(int j = 0; j < frame_idx.size(); j++){
                    int idx = frame_idx[j];
                    if (idx >= a1 && idx <= a2 )
                        meter++;
                }
                if (meter ==0)
                    KK--;

                meter = 0;

            }
        }

        global.codebook.init_k_for_kmeans += (KK+1);
        std::cout << "global.codebook.init_k_for_kmeans=" << global.codebook.init_k_for_kmeans<< std::endl;

    }
    


}


float calHammingDist(const cv::Mat& x1, const cv::Mat& x2){
    // x1 and x2 should be either binary or bipolar
    float dist = 0.0f;
    for(int i = 0; i < x1.cols; i++)
        dist += (x1.at<float>(i) == x2.at<float>(i)) ? 0 : 1.0f;

    return dist;
}


void jointsTemporalSmoothingKeyframeSelection() 
// Input: output_vec_KFrameBuffer
// Output: output_vec_KFramePool, output_vec_KFramePoolSkeletion (the skeletons after processing)
// Work: (1) temporal interpolation and (2) keyframe selection
{

    // (1) processing to handle oclusion and stablise tracking

    // if(global.n_frames_inBuffer < global.KEYFRAME_BUFFERSIZE){
    //     std::cout << "n_frames_in_buffer="<<global.n_frames_inBuffer<<std::endl;
    //     return nullptr;
    // }
    // int n_frames = global.KEYFRAME_BUFFERSIZE;
    // cv::Mat joints_o = cv::Mat::zeros(n_frames,global.num_dims,CV_32F);
    // cv::Mat joints = cv::Mat::zeros(n_frames,global.num_dims,CV_32F);
    // cv::Mat pose_pattern, extension, movement;
    // cv::Mat pose_pattern_bipolar;
    // cv::Mat memory_energy; // not necessary to be positive
    // for (int i = 0; i < n_frames; i++)
    //     global.output_vec_KFrameBuffer[i].copyTo(joints_o.row(i));
    
    // if (global.IS_DEBUG)
    //     std::cout << "[DEBUG] jointsTemporalSmoothingKeyframeSelection: begin to run" << std::endl;
    // // boost::timer::cpu_timer timer;

    // // temporalSpline(joints_o,joints);
    // joints = joints_o.clone();
    // if (global.IS_DEBUG)
    //     std::cout<< "break1" << std::endl;
    // calPosePattern (joints, pose_pattern, extension, movement);
    // if (global.IS_DEBUG)
    //     std::cout<< "break2" << std::endl;
    // calBipolarFeature(pose_pattern, pose_pattern_bipolar);
    // if (global.IS_DEBUG)
    //     std::cout<< "break3" << std::endl;
    // calBidirectionalAssociativeMemoryEnergy(pose_pattern_bipolar, extension, movement, memory_energy);
    // if (global.IS_DEBUG)
    //     std::cout<< "break4" << std::endl;
    // cv::normalize(memory_energy, memory_energy, 0,1,cv::NORM_MINMAX);
    // cv::normalize(movement, movement, 0,1, cv::NORM_MINMAX);
    // cv::GaussianBlur(memory_energy,memory_energy, cv::Size(1,25), 0, 12.5, cv::BORDER_REFLECT);
    // cv::GaussianBlur(movement,movement, cv::Size(1,25), 0, 12.5, cv::BORDER_REFLECT);
    // if (global.IS_DEBUG)
    //     std::cout<< "break5" << std::endl;



    // // if (global.IS_DEBUG){
    //     // std::cout<< "movement = " << movement.t() << std::endl;
    //     // std::cout<< "memory_energy = " << memory_energy.t() << std::endl;
    //     // std::cout<< "memory_energy = " << pose_pattern_bipolar.t() << std::endl;
    // // }
    

    // std::vector<int> frame_idx;    
    // getKeyFrameIdxFromBuffer(memory_energy, movement, frame_idx);
    // if (global.IS_DEBUG)
    //     std::cout<< "break6" << std::endl;

    /* --------------------result evaluation-------------------------------*/
    // std::cout <<" frame_idx.size()="<<frame_idx.size() <<std::endl;
    for (int i = 0; i < global.output_vec_KFrameBuffer.size(); i++) {
        global.output_vec_KFramePool.push_back(global.output_vec_KFrameBuffer[i]);
        // global.n_frames_inPool++;
        
        // global.output_vec_KFramePoolSkeleton.push_back(joints.row(frame_idx[i]+1));
    }
    if (global.IS_DEBUG)
        std::cout<< "break7" << std::endl;

    global.output_vec_KFrameBuffer.clear();


    // determine K for Kmeans. This module only works at the first working memory for initialisation
    // However, the value of K is not equivalent to the temporal segmentation. Therefore, we create a fully connected graph, where variables are bipolar and distance a hamming
    // if (global.KEYFRAME_POOL_FILLTIMES==0){
        
    //     // extract maxima in the memory+movement curve, so as to segment the video.
    //     cv::Mat _W (cv::Size(pose_pattern_bipolar.rows, pose_pattern_bipolar.rows), CV_32F, cv::Scalar(0));
    //     cv::Mat _D;
    //     cv::Mat _L;
    //     cv::Mat _L_eigenval;
    //     cv::Mat _aux_mat;
    //     int n_clusters;
        
    //     for (int i = 0; i < _W.rows; i++) {
    //         for (int j = i; j < _W.cols; j++){
    //             cv::Mat _jts1 = pose_pattern_bipolar.row(i);
    //             cv::Mat _jts2 = pose_pattern_bipolar.row(j);
    //             _W.at<float> (i,j) = calHammingDist(_jts1, _jts2); // 
    //         }
    //     } 

    //     _W = _W + _W.t() - cv::Mat::diag(_W.diag());
    //     // std::cout << "[DEBUG] _W = " << _W <<std::endl;
    //     cv::reduce(_W, _D, 0, CV_REDUCE_SUM, CV_32F);

    //     _L = cv::Mat::diag(_D)-_W; // unnormalized graph laplacian
    //     // _L =  cv::Mat::eye(_W.size(), CV_32F) - cv::Mat::diag(_D).inv() * _W; // normalized graph laplacian. See the literature

    //     std::cout << "L = " << _L << std::endl;
    //     cv::eigen(_L, _L_eigenval);
    //     cv::sort(_L_eigenval, _L_eigenval, CV_SORT_EVERY_COLUMN+CV_SORT_ASCENDING);
    //     cv::sortIdx( _L_eigenval(cv::Range(1,_W.rows),cv::Range::all())-_L_eigenval(cv::Range(0,_W.rows-1),cv::Range::all()), 
    //         _aux_mat, CV_SORT_EVERY_COLUMN+CV_SORT_DESCENDING);
    //     n_clusters = _aux_mat.at<int>(0)+1;    

    //     global.codebook.init_k_for_kmeans +=  _aux_mat.at<int>(0)+1;
    //     std::cout << "global.codebook.init_k_for_kmeans=" << global.codebook.init_k_for_kmeans<< std::endl;

    // }



}







void initCodebook ( ) {
    // std::cout << "[Codebook Initialization] starts.." << std::endl;

    if (global.VERBOSE)
        std::cout << "[CODEBOOK INFO INIT] compute from data..." << std::endl;
    
    // 1. data standardization
    cv::Mat samples_o = cv::Mat::zeros(global.output_vec_KFramePool.size(), global.num_dims, CV_32F); // the detection score is included
    cv::Mat samples;
    for (int i = 0; i < global.output_vec_KFramePool.size(); i++) 
       global.output_vec_KFramePool[i].copyTo(samples_o.row(i));
       
    // std::cout << "[DEBUG] samples matrix size= " << samples_o.size() << std::endl;
    sampleStandardization(samples_o, samples, global.USE_STANDARDIDATION);

    if (global.KEYFRAME_BUFFERSIZE == 1)
        global.codebook.init_k_for_kmeans = 1;
    else
        global.codebook.init_k_for_kmeans = calClusterNumbers(samples);


    global.codebook.n_clusters = global.codebook.init_k_for_kmeans;
    
    if (global.VERBOSE){
        std::cout << "[Codebook Initialization] #clusters=" << global.codebook.n_clusters << std::endl;
        std::cout<< "[Codebook Initialization] #samples=" << global.output_vec_KFramePool.size() <<std::endl;
    }
    // 3. according to n_clusters, we perform kmeans to specify variables in global.codebook
    
    // cv::Mat samples;
    cv::Mat labels, cluster_locs;
    // std::cout << "[DEBUG] global.codebook.init_k_for_kmeans="<<global.codebook.init_k_for_kmeans <<std::endl;
    if (global.codebook.init_k_for_kmeans > 1){
        cv::kmeans(samples, global.codebook.init_k_for_kmeans, labels, cv::TermCriteria(cv::TermCriteria::COUNT+cv::TermCriteria::EPS, 1000,1e-6),
               1, cv::KMEANS_PP_CENTERS, cluster_locs);
    }
    else{
        labels = cv::Mat::zeros(samples.rows, 1, CV_32F);
        cv::reduce(samples, cluster_locs, 0, CV_REDUCE_AVG);
    }
    labels.convertTo(labels, CV_32F);
    
    // std::ofstream fout ("samples.txt");
    // fout << samples;
    // std::cout<< "[DEBUG] labels =" << labels <<std::endl;
    global.codebook.init_sample_labels = labels.clone();
    global.codebook.cluster_locs = cluster_locs.clone();
    for(int i = 0; i < samples.rows; i++)
        global.codebook.sample_labels.push_back(labels.at<float>(i));

    for(int i = 0; i < global.codebook.init_k_for_kmeans; i++)
        global.codebook.cluster_labels.push_back(i);

    if (global.VERBOSE){
        std::cout <<"--compute cluster radius" << std::endl;
    }
    calClusterStd(samples);
    // calClusterImgs();
    // intialize the transition map to the uniform distribution( we save all the elements to 1)

    if (global.VERBOSE){
        std::cout <<"--compute cluster radius finish" << std::endl;
    }
    cv::Mat transition_map_init = cv::Mat::ones(global.codebook.init_k_for_kmeans, global.codebook.init_k_for_kmeans, CV_32F);

    // transition_map_init = transition_map_init / cv::sum(transition_map)[0];
    global.codebook.transition_map = cv::SparseMat(transition_map_init);
    global.codebook.index_cluster_pre = labels.at<int>(samples.rows-1);
    global.codebook.prob_prediction = cv::Mat::ones(global.codebook.init_k_for_kmeans,1, CV_32F);        
   
    if (global.VERBOSE){
        std::cout << "[CODEBOOK INFO INIT] codebook.n_clusters = " << global.codebook.n_clusters << std::endl;
        std::cout << "[CODEBOOK INFO INIT] codebook.cluster_labels = " << global.codebook.cluster_labels.t() << std::endl;
        std::cout << "[CODEBOOK INFO INIT] codebook.cluster_locs = " << global.codebook.cluster_locs << std::endl;
        std::cout << "[CODEBOOK INFO INIT] codebook.cluster_std = " << global.codebook.cluster_std<< std::endl;
        cv::Mat n_samples_per_cluster = cv::Mat(global.codebook.n_samples_per_cluster);
        std::cout << "[CODEBOOK INFO INIT] codebook.n_samples_per_cluster = " << n_samples_per_cluster.t()<< std::endl;
        std::cout << "[CODEBOOK INFO INIT] codebook.prob_prediction = "<< global.codebook.prob_prediction.t() << std::endl;
    }
}





void calDistSampleToCluster(cv::Mat& phi, cv::Mat& cluster_locs, cv::Mat& dist) {
    
    for(int j = 0; j < cluster_locs.rows; j++) {
        float _dist = 0.0f;

        for (int i = 0; i < phi.cols; i++){
            float x1 = phi.at<float> (i);
            float x2 = cluster_locs.at<float> (j,i);
            
            // float ss = (x1 == 0.0f || x2 == 0.0f || y1 == 0.0f || y2 == 0.0f)? 0.0f : 1.0f;
            _dist += (x1-x2)*(x1-x2) ;
        }
        dist.push_back(std::sqrt(_dist));
    }
}


// p(C_j | x) = \int_{C_i} p(C_j | C_i)p(C_i | x) 
void activityPredict(cv::Mat& x, cv::Mat& prob_prediction) {


    // int n_dims = global.codebook.cluster_locs.cols;
    int n_clusters = global.codebook.n_clusters;
    cv::Mat prob_transition = cv::Mat::zeros(n_clusters,n_clusters, CV_32F);
    cv::Mat prob_observation(1,n_clusters, CV_32F);

    // calculate p(C_i | x)
    for(int i = 0; i < n_clusters; i++){
        // cv::Mat _exp_term;

        cv::Mat mu = global.codebook.cluster_locs.row(i);
        // cv::Mat sigma2 = (global.codebook.cluster_std.row(i)).mul(global.codebook.cluster_std.row(i));
        // std::cout << "[DEBUG] activityPredict(): mu = " << mu << std::endl;
        // std::cout << "[DEBUG] activityPredict(): sigma2 = " << mu << std::endl;

        // cv::exp(-((x-mu).mul(x-mu)).mul( 1.0f/(2.0f*sigma2) ), _exp_term);
        // cv::Mat _Z;
        // cv::sqrt(1.0f/(2*3.1415926*sigma2), _Z);
        // cv::Mat _g, _g_one;
        // _g = _exp_term.mul(_Z);
        // cv::log(_g,_g);
        // cv::reduce(_g, _g_one, 0, CV_REDUCE_SUM);
        // cv::exp(_g_one,_g_one);
        prob_observation.at<float>(i) = cv::exp(-0.001f*cv::norm(x,mu,cv::NORM_L2));
        // std::cout << "[DEBUG] activityPredict(): cv::norm(x,mu,cv::NORM_L2) = " << cv::norm(x,mu,cv::NORM_L2) <<std::endl;
    }
    // std::cout << "[DEBUG] activityPredict(): prob_observation_before=" << prob_observation << std::endl;

    cv::normalize(prob_observation, prob_observation,1,0,cv::NORM_L1);
    // std::cout << "[DEBUG] activityPredict(): prob_observation_after=" << prob_observation << std::endl;
    
    // std::cout << "[DEBUG] activityPredict(): current prob_observation = " << prob_observation << std::endl;
    // calculate p(C_j | C_i), convert sparseMat to Mat
    cv::SparseMat _p_tran = global.codebook.transition_map.clone();
    // int dims = _p_tran.dims();
    
    // std::cout << "[DEBUG] activityPredict(): create sparsemat iterator" << std::endl;
    cv::SparseMatConstIterator_<float> 
                it = _p_tran.begin<float>(),
                it_end = _p_tran.end<float>();

    // std::cout << "[DEBUG] activityPredict(): assign sparse to dense mat" << std::endl;
    for(; it != it_end; ++it) {

        // const cv::SparseMat::Node* n = it.node();
        prob_transition.at<float>(it.node()->idx[0],it.node()->idx[1]) =  _p_tran.ref<float>(it.node()->idx);
        // std::cout << "[DEBUG] activityPredict(): probt_ransition.at<float>(index[0],index[1]) = " << prob_transition.at<float>(it.node()->idx) << std::endl;    
    }
    // std::cout << "[DEBUG] activityPredict(): current prob_transition = " << prob_transition << std::endl;
    prob_prediction = prob_transition.t()*prob_observation.t(); // column vector
    // global.codebook.transition_map.clear();
    // global.codebook.transition_map = cv::SparseMat(prob_transition);

    cv::normalize(prob_prediction, prob_prediction,1,0,cv::NORM_L1);

}




void updateCodebook() {
    for(int i = 0; i < global.output_vec_KFramePool.size(); i++) {
        cv::Mat _dist, _joints, _joints_t, _joints_repmat;
        double _dist_min;
        int _dist_min_idx; 
        
        std::vector<float> _cluster_radius;
        _joints = global.output_vec_KFramePool[i];
        // std::cout << _joints << std::endl;
        // _joints = global.output_vec_KFramePoolSkeleton[i];
        
        // _joints_t = _joints.t();
        // cv::copyTo(_joints.t(), _joints_t);
        // cv::transpose(_joints, _joints_t);
        // cv::divide((_joints_t - global.codebook.init_standardization_mean), global.codebook.init_standardization_std), _joints_t);
        _joints_t = (_joints - global.codebook.init_standardization_mean) / global.codebook.init_standardization_std;

        // _joints_t = _joints.t();
        // cv::repeat(_joints_t,  global.codebook.n_clusters, 1, _joints_repmat);
        // // std::cout << _joints_repmat << std::endl;
        // // cv::waitKey(0);
        // cv::reduce((global.codebook.cluster_locs - _joints_repmat).mul(global.codebook.cluster_locs - _joints_repmat),
        //     _dist2, 1, CV_REDUCE_SUM);
        // std::cout << "[DEBUG] _joints_t = " << _joints_t << std::endl;
        calDistSampleToCluster(_joints_t, global.codebook.cluster_locs, _dist);
        cv::minMaxIdx(_dist, &_dist_min, NULL, &_dist_min_idx);
        calClusterRadiusfromStd(_cluster_radius);
        float _radius_max = *std::max_element(_cluster_radius.begin(), _cluster_radius.end() );
        // std::cout<< "[DEBUG] updateCodebook(): _dist_min = " << _dist_min << ", _radius_max = " << _radius_max << std::endl;
        // std::cout<< "[DEBUG] updateCodebook(): _dist = " << _dist.t() <<std::endl;
        // std::cout<< "[DEBUG] updateCodebook(): _dist_min_idx = " << _dist_min_idx <<std::endl;


        // check + global.RADIUS_OFFSET transitional abnormality and update activity prediction probility
        // cv::sortIdx(global.codebook.prob_prediction, idx_prob_prediction, CV_SORT_EVERY_COLUMN+CV_SORT_DESCENDING);
        // if (global.codebook.prob_prediction.at<float>(_dist_min_idx) < 0.5*1.0f/global.codebook.n_clusters)
            // std::cout<< "[CODEBOOK INFO] [ABNORMAL!!!] dynamical abnormal occurs!" << std::endl;
         

        // std::cout << "--------------------------------------------------------------------------------" << std::endl;
        // if(global.IS_DEBUG){
            // std::cout << "[DEBUG] _dist_min = " << _dist_min << "  [DEBUG] _radius_max = " << _radius_max << std::endl;
        // }
        if (global.VERBOSE)
            std::cout << "[DEBUG] _dist_min = " << _dist_min << "  [DEBUG] _radius_max = " << _radius_max << std::endl;
        
        if ( _dist_min > _radius_max) { // create a new cluster with new location, std=0 and image
            if (global.VERBOSE)
                std::cout<< "[CODEBOOK INFO] updateCodebook(): creating new cluster";
            // std::cout<< "[CODEBOOK INFO] [ABNORMAL!!!] abnormal appearance occurs!" << std::endl;
            global.codebook.cluster_labels.push_back(global.codebook.n_clusters);
            global.codebook.sample_labels.push_back((float)global.codebook.n_clusters);
            global.codebook.n_clusters += 1;
            // std::cout<< "[DEBUG] updateCodebook(): creating new cluster, mean" << std::endl;
            global.codebook.cluster_locs.push_back(_joints_t);
            // std::cout<< "[DEBUG] updateCodebook(): creating new cluster, std" << std::endl;
            cv::Mat zz = cv::Mat::zeros(1,_joints_t.cols, CV_32F) + global.codebook.init_cluster_std; // avoid numerical singularity
            global.codebook.cluster_std.push_back(zz);
            // std::cout<< "[DEBUG] updateCodebook(): creating new cluster, frame" << std::endl;
            // global.codebook.cluster_imgs.push_back(_frame);
            // std::cout<< "[DEBUG] updateCodebook(): creating new cluster, #asmples" << std::endl;
            global.codebook.n_samples_per_cluster.push_back(1);

            // update transition map
            // std::cout<< "[DEBUG] updateCodebook(): creating new cluster, create new transition node" << std::endl;

            int transition_map_idx[global.codebook.transition_map_dims];
            transition_map_idx[0] = global.codebook.index_cluster_pre;
            transition_map_idx[1] = global.codebook.n_clusters-1;
            global.codebook.transition_map.newNode(transition_map_idx, 0);
            global.codebook.transition_map.ref<float>(transition_map_idx) = 1.0f;
            // std::cout<< "[DEBUG] updateCodebook(): creating new cluster, create new transition node end" << std::endl;

            // float* _tran_ptr = reinterpret_cast<float*>(global.codebook.transition_map.ptr(transition_map_idx[0],transition_map_idx[1], true));
            // *_tran_ptr = 1.0f;
            // global.codebook.transition_map.addref();
            // global.codebook.transition_map.ref<float>(transition_map_idx[0], transition_map_idx[1]) = 1.0f;
            
            global.codebook.index_cluster_pre = transition_map_idx[1];
            // std::cout<< "[DEBUG] updateCodebook(): creating new cluster, create new transition node end 2" << std::endl;


        }
        else {  // update the cluster and the cluster image, according to _dist2_min_idx
            
            // std::cout << "[CODEBOOK INFO CURRENT OBSERVATION] O(C_t) = " << _dist_min_idx << std::endl;
            // std::cout << "[CODEBOOK INFO PREVIOUS PREDICTION] P(" << _dist_min_idx << "| x_{t-1}) = " 
                      // << global.codebook.prob_prediction.at<float>(_dist_min_idx) << std::endl;
            // if (global.codebook.prob_prediction.at<float>(_dist_min_idx) < 0.05)
            //     std::cout<< "[CODEBOOK INFO] [ABNORMAL!!!] abnormal movement occurs!" << std::endl;
            if (global.VERBOSE)
                std::cout<< "[CODEBOOK INFO] updateCodebook(): updating current cluster " << _dist_min_idx << std::endl;
            global.codebook.sample_labels.push_back((float)_dist_min_idx);
            int _n_samplei = global.codebook.n_samples_per_cluster[_dist_min_idx];
            cv::Mat _mean = global.codebook.cluster_locs.row(_dist_min_idx).clone();
            cv::Mat _std = global.codebook.cluster_std.row(_dist_min_idx);
            // updating mean
            // std::cout<< "[DEBUG] updateCodebook(): updating mean " << std::endl;

            global.codebook.cluster_locs.row(_dist_min_idx) = _mean + (_joints_t-_mean ) / (float) (1+_n_samplei);

            // updating std
            // std::cout<< "[DEBUG] updateCodebook(): updating std " << std::endl;
            cv::Mat _var;
            cv::pow(_std, 2, _var);
            _var = ((_joints_t-_mean).mul(_joints_t-global.codebook.cluster_locs.row(_dist_min_idx)) + _n_samplei*_var)/(_n_samplei+1);
            cv::sqrt(_var, global.codebook.cluster_std.row(_dist_min_idx));

            // updating transition_map
            int transition_map_idx[global.codebook.transition_map_dims];
            transition_map_idx[0] = global.codebook.index_cluster_pre;
            transition_map_idx[1] = _dist_min_idx;
            global.codebook.transition_map.ref<float>(transition_map_idx)+=1.0f;

            // updating cluster image
            // std::cout<< "[DEBUG] updateCodebook(): updating img "<< std::endl;
            // cv::Mat _cluster_img_joint;
            // cvtFrameToJoints(global.codebook.cluster_imgs[_dist_min_idx], _cluster_img_joint);
            // // _cluster_img_joint = global.
            // // std::cout<< "[DEBUG] updateCodebook(): _cluster_img_joint=" << _cluster_img_joint<< std::endl;
            // float _dist_img2center = cv::norm(_cluster_img_joint.t(), _mean);
            // if (_dist_img2center > (_dist_min) )
            //     global.codebook.cluster_imgs[_dist_min_idx] = _frame;

            // update sample numbers
            global.codebook.n_samples_per_cluster[_dist_min_idx] = _n_samplei + 1;
            // std::cout << "[CODEBOOK INFO] updateCodebook(): current num_samples = " << global.codebook.n_samples_per_cluster[_dist_min_idx]<< std::endl;    // std::cout << "imshowing starting 2" << std::endl;
            global.codebook.index_cluster_pre = _dist_min_idx;
        }

        // std::cout<< "[DEBUG] updateCodebook(): creating new cluster, perform prediction" << std::endl;
        // activityPredict(_joints_t, global.codebook.prob_prediction);
        // std::cout<< "[DEBUG] updateCodebook(): creating new cluster, perform prediction end" << std::endl;
    }
    // std::cout << "sample_labels = " << global.codebook.sample_labels.rowRange(global.codebook.sample_labels.rows-global.output_vec_KFramePool.size(), global.codebook.sample_labels.rows).t() << std::endl;
    // std::cout << "global.codebook.n_clusters = " << global.codebook.n_clusters << std::endl;
    // std::cout << "global.codebook.cluster_locs = " << global.codebook.cluster_locs << std::endl;
    // std::cout << "global.codebook.cluster_std = " << global.codebook.cluster_locs << std::endl;
    // std::cout << "global.codebook.n_samples_per_cluster = " << global.codebook.n_samples_per_cluster << std::endl;
}





void updateCodebookTemporalRegularity() {

    int n_clusters = global.codebook.n_clusters;
    std::vector<int> n_samples_per_cluster = global.codebook.n_samples_per_cluster;
    cv::Mat cluster_labels = global.codebook.cluster_labels.clone();
    cv::Mat cluster_locs = global.codebook.cluster_locs.clone();
    cv::Mat cluster_std = global.codebook.cluster_std.clone();
    // cv::Mat sample_labels = global.codebook.sample_labels.clone();
    cv::Mat sample_labels;

    // the first round, assign new labels to new samples.
    for(int i = 0; i < global.output_vec_KFramePool.size(); i++) {
        cv::Mat _dist, _joints, _joints_t, _joints_repmat;
        double _dist_min;
        int _dist_min_idx; 
        
        std::vector<float> _cluster_radius;
        _joints = global.output_vec_KFramePool[i];
        _joints_t = (_joints - global.codebook.init_standardization_mean) / global.codebook.init_standardization_std;

        calDistSampleToCluster(_joints_t, cluster_locs, _dist);
        cv::minMaxIdx(_dist, &_dist_min, NULL, &_dist_min_idx);
        calClusterRadiusfromStd(_cluster_radius);
        float _radius_max = *std::max_element(_cluster_radius.begin(), _cluster_radius.end() );
        

        
        if ( _dist_min > _radius_max) { // create a new cluster with new location, std=0 and image
            cluster_labels.push_back(n_clusters);
            sample_labels.push_back((float)n_clusters);
            n_clusters += 1;
            cluster_locs.push_back(_joints_t);
            cv::Mat zz = cv::Mat::zeros(1,_joints_t.cols, CV_32F) + global.codebook.init_cluster_std; // avoid numerical singularity
            cluster_std.push_back(zz);
            n_samples_per_cluster.push_back(1);

        }
        else {  

            int _n_samplei = n_samples_per_cluster[_dist_min_idx];
            sample_labels.push_back((float)_dist_min_idx);
            cv::Mat _mean = cluster_locs.row(_dist_min_idx);
            cv::Mat _std = cluster_std.row(_dist_min_idx);

            cluster_locs.row(_dist_min_idx) = _mean + (_joints_t-_mean ) / (float) (1+_n_samplei);

            cv::Mat _var;
            cv::pow(_std, 2, _var);
            _var = ((_joints_t-_mean).mul(_joints_t-_mean) + _n_samplei*_var)/(_n_samplei+1);
            cv::sqrt(_var, cluster_std.row(_dist_min_idx));

            n_samples_per_cluster[_dist_min_idx] = _n_samplei + 1;            
        }
    }


    // temporal regularity
    cv::Mat sample_labels_t, sample_labels_sm;
    // int N = sample_labels.rows;
    // sample_labels_t = sample_labels.rowRange(N-global.output_vec_KFramePool.size(), N);
    // sample_labels_t = sample_labels.t();
    temporalSpline1D(sample_labels, sample_labels_sm);
    // std::cout << "sample_labels = " << sample_labels.t() << std::endl;
    // std::cout << "sample_labels_sm = " << sample_labels_sm.t() << std::endl;


    // update codebook depending on temporal smoothed label
    for(int i = 0; i < global.output_vec_KFramePool.size(); i++) {
    
        int label = (int)roundf(std::abs(sample_labels_sm.at<float>(i)));
        cv::Mat feature = global.output_vec_KFramePool[i];
        global.codebook.sample_labels.push_back((float)label);
        // if this label is new
        if (label >= global.codebook.n_clusters){
            global.codebook.n_clusters ++;
            global.codebook.cluster_labels.push_back(label);
            global.codebook.cluster_locs.push_back(feature);

            cv::Mat zz = cv::Mat::zeros(1,feature.cols, CV_32F) + global.codebook.init_cluster_std; // avoid numerical singularity
            global.codebook.cluster_std.push_back(zz);

            global.codebook.n_samples_per_cluster.push_back(1);
            
            // update transition map
            int transition_map_idx[global.codebook.transition_map_dims];
            transition_map_idx[0] = global.codebook.index_cluster_pre;
            transition_map_idx[1] = global.codebook.n_clusters-1;
            global.codebook.transition_map.newNode(transition_map_idx, 0);
            global.codebook.transition_map.ref<float>(transition_map_idx) = 1.0f;
            // float* _tran_ptr = reinterpret_cast<float*>(global.codebook.transition_map.ptr(transition_map_idx[0],transition_map_idx[1], true));
            // *_tran_ptr = 1.0f;
            // global.codebook.transition_map.addref();
            // global.codebook.transition_map.ref<float>(transition_map_idx[0], transition_map_idx[1]) = 1.0f;
            
            global.codebook.index_cluster_pre = transition_map_idx[1];

        }
        else{
        // else, we update the corresponding cluster and the transition matrix.
            int _n_samplei = global.codebook.n_samples_per_cluster[label];
            
            cv::Mat _mean = global.codebook.cluster_locs.row(label);
            cv::Mat _std = global.codebook.cluster_std.row(label);
            // updating mean
            // std::cout<< "[DEBUG] updateCodebook(): updating mean " << std::endl;

            global.codebook.cluster_locs.row(label) = _mean + (feature-_mean ) / (float) (1+_n_samplei);

            // updating std
            // std::cout<< "[DEBUG] updateCodebook(): updating std " << std::endl;
            cv::Mat _var;
            cv::pow(_std, 2, _var);
            _var = ((feature-_mean).mul(feature-_mean) + _n_samplei*_var)/(_n_samplei+1);
            cv::sqrt(_var, global.codebook.cluster_std.row(label));

            // updating transition_map
            int transition_map_idx[global.codebook.transition_map_dims];
            transition_map_idx[0] = global.codebook.index_cluster_pre;
            transition_map_idx[1] = label;
            global.codebook.transition_map.ref<float>(transition_map_idx)+=1.0f;
            global.codebook.index_cluster_pre = label;
            
            // update sample numbers
            global.codebook.n_samples_per_cluster[label] = _n_samplei + 1;
        }
        activityPredict(feature, global.codebook.prob_prediction);

    }

    // std::cout << "global.codebook.n_clusters = " << global.codebook.n_clusters << std::endl;
    // std::cout << "global.codebook.cluster_locs = " << global.codebook.cluster_locs << std::endl;
    // std::cout << "global.codebook.cluster_std = " << global.codebook.cluster_locs << std::endl;
    // std::cout << "global.codebook.n_samples_per_cluster = " << global.codebook.n_samples_per_cluster << std::endl;
    
    
    sample_labels.release();
    sample_labels_t.release();
    sample_labels_sm.release();

}








// in this function, keyframe candidate selection and incremental clustering run successively
void runIncrementalClustering() {

    // 1. keyframe candidate selection, the results are stored in global.output_vec_KFramePool
    // NOTICE: in the current implementation, this function is disabled.
        
    if (global.IS_DEBUG)
        std::cout << "[DEBUG] extract frames from working memory" << std::endl;
    jointsTemporalSmoothingKeyframeSelection();


    // 2. incremental learning
    
// if (global.output_vec_KFramePool.size() >= global.MIN_KEYFRAME_POOLSIZE){
    if (global.VERBOSE){
        std::cout << "[CODEBOOK INFO] global.KEYFRAME_POOL_FILLTIMES=" << global.KEYFRAME_POOL_FILLTIMES << std::endl;
        std::cout << "[CODEBOOK INFO] global.KEYFRAME_POOL_SIZE=" << global.output_vec_KFramePool.size() << std::endl;
    }
    
    if (global.KEYFRAME_POOL_FILLTIMES ==0) { // initialization using kmeansr
        if (global.VERBOSE)
            std::cout << "[CODEBOOK INFO] runIncrementalClglobal.n_featurestering(): Codebook Initialization--------------------BEGIN" << std::endl;
        initCodebook( );
        if (global.VERBOSE)
            std::cout << "[CODEBOOK INFO] runIncrementalClustering(): Codebook Initialization--------------------END" << std::endl;

    }
    else {
        if (global.VERBOSE)
            std::cout << "[CODEBOOK INFO] runIncrementalClustering(): Codebook Updating--------------------BEGIN" << std::endl;
        if(!global.USE_TEMPORAL_REG)
            updateCodebook( );
        else
            updateCodebookTemporalRegularity();
        if (global.VERBOSE)
            std::cout << "[CODEBOOK INFO] runIncrementalClustering(): Codebook Updating--------------------END" << std::endl;
    }
    
    global.KEYFRAME_POOL_FILLTIMES+=1;
    global.output_vec_KFramePool.clear();
    // }

    

}





void assignLabelsToSamples(const cv::Mat& samples, cv::Mat& labels){
    double _dist_min;
    int _dist_min_idx;

    for (int i = 0; i < global.num_samples; i++){
        cv::Mat phi = samples.row(i);
        cv::Mat dist;

        calDistSampleToCluster(phi, global.codebook.cluster_locs, dist);
        cv::minMaxIdx(dist, &_dist_min, NULL, &_dist_min_idx);
        labels.push_back(_dist_min_idx);
    }
}


// Input: all the pose patterns. Rows - samples, cols - features
// Output: codebook, segment boudnaries.
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray *prhs[])
{

    if (nrhs != 6)
        mexErrMsgIdAndTxt("Three arguments should be parsed", "usage: fun(samples, ShortMemorySize, errorEstimate, USE_TEMPORAL_REG, VERBOSE, TVL1_weights)");

    // if (!mxIsDouble(prhs[0]))
    //     mexErrMsgIdAndTxt("The first element is not double", "usage: fun(samples, ShortMemorySize, errorEstimate, USE_TEMPORAL_REG)");

    // if (!mxIsDouble(prhs[1]))
    //     mexErrMsgIdAndTxt("The second element is not double", "usage: fun(samples, ShortMemorySize, errorEstimate, USE_TEMPORAL_REG)");


    // parse arguments
    global.KEYFRAME_BUFFERSIZE = (int)*mxGetPr(prhs[1]);
    int offset = global.KEYFRAME_BUFFERSIZE;
    float error_estimate = (float)*mxGetPr(prhs[2]);
    global.USE_TEMPORAL_REG = (int)*mxGetPr(prhs[3]);
    if (global.KEYFRAME_BUFFERSIZE != 1)
        global.USE_TEMPORAL_REG = (int)*mxGetPr(prhs[3]);
    else
        global.USE_TEMPORAL_REG = 0;
    global.VERBOSE = (int)*mxGetPr(prhs[4]);
    global.temporal_reg_alpha = (float)*mxGetPr(prhs[5]);

    cv::Mat all_samples;
    ocvMxArrayToMat_double(prhs[0], all_samples);
    all_samples.convertTo(all_samples, CV_32F);
    
    global.num_samples = all_samples.rows;
    global.num_dims = all_samples.cols;
    global.KEYFRAME_POOL_FILLTIMES = 0;
    global.RADIUS_OFFSET = error_estimate*std::sqrt(global.num_dims);
    global.output_vec_KFrameBuffer.clear();
    global.output_vec_KFramePool.clear();
    global.codebook.cluster_labels.release();
    global.codebook.sample_labels.release();
    // step2: run the incremental learning mechanism
    for(int i = 0; i < global.num_samples; i+=offset){
        
        int upper_idx_bound = std::min(i+global.KEYFRAME_BUFFERSIZE, global.num_samples);
        // std::cout << upper_idx_bound << std::endl;
        // std::cout << global.num_samples << std::endl;
        // std::cout << i << std::endl;
        cv::Mat samples_for_buffer = all_samples.rowRange(i, upper_idx_bound);
        for(int j = 0; j < samples_for_buffer.rows; j++)
            global.output_vec_KFrameBuffer.push_back(samples_for_buffer.row(j));
        runIncrementalClustering();
    }

    // showCodebookInfo();
    // assign the labels to all samples
    // cv::Mat labels, labels_tr;
    // // assignLabelsToSamples(all_samples, labels);
    // labels.convertTo(labels, CV_32F);
    // std::cout << "perform temporal interjpoatioin" << std::endl;
    // std::cout << labels << std::endl;
    // temporalSpline1D(labels,labels_tr);

    // std::cout << "perform temporal interjpoatioin. Finish" << std::endl;
    // labels_tr.convertTo(labels_tr, CV_64F);
    // labels.convertTo(labels, CV_64F);
    cv::Mat labels;
    cv::Mat cluster_locs;
    global.codebook.sample_labels.convertTo(labels, CV_64F);
    global.codebook.cluster_locs.convertTo(cluster_locs, CV_64F);
    if (global.VERBOSE)
        std::cout << "before write to mxArray" << std::endl;
    // plhs[1] = ocvMxArrayFromMat_double(labels_tr);
    

    // cv::Mat prob_transition = cv::Mat::zeros(global.codebook.n_clusters,global.codebook.n_clusters,CV_32F);
    // cv::SparseMatConstIterator_<float> 
    //             it = global.codebook.transition_map.begin<float>(),
    //             it_end = global.codebook.transition_map.end<float>();

    // for(; it != it_end; ++it) {

    //     const cv::SparseMat::Node* n = it.node();
    //     const int* index = n->idx;
    //     prob_transition.at<float>(index[0],index[1]) =  global.codebook.transition_map.ref<float>(index);
    //     // std::cout << "[DEBUG] activityPredict(): prob_transition.at<float>(index[0],index[1]) = " << prob_transition.at<float>(index[0],index[1]) << std::endl;    
    // }
    // prob_transition.convertTo(prob_transition, CV_64F);

    plhs[0] = ocvMxArrayFromMat_double(labels);
    plhs[1] = ocvMxArrayFromMat_double(cluster_locs);
    // plhs[2] = ocvMxArrayFromMat_double(prob_transition);

}

