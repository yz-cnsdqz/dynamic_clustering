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


      
// global queues for I/O
struct Global {

    int num_samples = 0;
    int num_dims = 0;
    int KEYFRAME_BUFFERSIZE = 0;       // BUFFER size. Reaching this size triggers [ALGO: keyframe selection] and then release the buffer.
    
    int KEYFRAME_POOL_FILLTIMES = 0;    // when the keyframe pool is full, then perform clustering and then released, this number +1.
    static const bool USE_STANDARDIDATION = false;
    
    int VERBOSE = 0;
    // static const bool SAVE_IMAGE = true;
    // static const bool SAVE_CODEBOOK=true;
    std::vector<cv::Mat> output_vec_KFrameBuffer; // a set of frames are stored within the buffer ( not always from the beginning)
    std::vector<cv::Mat> output_vec_KFramePool;   // the set of keyframe candidiates. The sorted result
    static const bool IS_DEBUG = false;
    float delta = 0.0f; // compensate floating error
    float theta = 1.0f;
    struct Codebook {
        int n_clusters;
        cv::Mat init_standardization_mean;
        cv::Mat init_standardization_std;
        cv::Mat init_sample_labels; 
        cv::Mat sample_labels;
        
        float init_cluster_std = 0.0f;
        int init_k_for_kmeans = 0;
        std::vector<int> n_samples_per_cluster;
        
        cv::Mat cluster_labels;
        
        cv::Mat cluster_locs; 
        cv::Mat cluster_std; 
        cv::Mat cluster_Ex2; 
        

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
        cv::Mat _var = _std.mul(_std);
    
        float  _var_sum = cv::sum(_var)[0];
        float r_i = (global.theta * _var_sum>=global.delta)? global.theta * _var_sum : global.delta;
        radius.push_back(r_i);
    }
}

void calClusterStd(cv::Mat& samples) {
    
    global.codebook.n_samples_per_cluster = std::vector<int> (global.codebook.n_clusters, 0); 
    global.codebook.cluster_std = cv::Mat(global.codebook.cluster_locs.size(), CV_32F);
    global.codebook.cluster_Ex2 = cv::Mat(global.codebook.cluster_locs.size(), CV_32F);
    
    for (int i = 0; i < global.codebook.n_clusters; i++){
        cv::Mat _samplei;       

        for (int j = 0; j < samples.rows; j++) {
            
            if (global.codebook.init_sample_labels.at<float>(j)==(float)i){
                _samplei.push_back(samples.row(j));   
                global.codebook.n_samples_per_cluster[i] +=1;
            }

        }
        
        int n_samplei = global.codebook.n_samples_per_cluster[i];
        
        cv::Mat _mean = global.codebook.cluster_locs.row(i);
        cv::Mat _mean_repmat = cv::repeat(_mean, n_samplei, 1);
        
        cv::Mat _std, _ex2;
        cv::reduce ((_samplei-_mean_repmat).mul((_samplei-_mean_repmat))/n_samplei, _std, 0, CV_REDUCE_SUM);
        cv::sqrt(_std,_std);
        _std.copyTo(global.codebook.cluster_std.row(i));
        

        cv::reduce ((_samplei).mul((_samplei))/n_samplei, _ex2, 0, CV_REDUCE_SUM);
        _ex2.copyTo(global.codebook.cluster_Ex2.row(i));
        
    }
}


void sampleStandardization(cv::Mat& samples, cv::Mat& samples_dst, bool flag) {
    
    if (flag){
    
        cv::Mat _mean, _mean_repmat;
        cv::Mat _std, _std_repmat;
        cv::reduce(samples, _mean, 0, CV_REDUCE_AVG);
    
        cv::repeat(_mean, samples.rows, 1, _mean_repmat);
        cv::reduce((samples-_mean_repmat).mul(samples-_mean_repmat), _std, 0, CV_REDUCE_AVG);
        cv::sqrt(_std, _std);
        cv::repeat(_std, samples.rows, 1, _std_repmat);
        global.codebook.init_standardization_mean = _mean.clone();
        global.codebook.init_standardization_std = _std.clone();    
    
        cv::divide(samples-_mean_repmat, _std_repmat, samples_dst);
    }
    else{
        samples_dst = samples.clone();
        global.codebook.init_standardization_mean = cv::Mat::zeros(1,global.num_dims, CV_32F);
        global.codebook.init_standardization_std = cv::Mat::ones(1,global.num_dims,CV_32F);
    }
    
}


void showCodebookInfo( ) {
    std::cout << "------------------------Final Result ---------------------------" <<std::endl;
    std::cout << "[CODEBOOK INFO] codebook.n_clusters = " << global.codebook.n_clusters << std::endl;
    
    // cv::namedWindow("Codebook Image");
    for (int k = 0; k < global.codebook.cluster_locs.rows; k++) {
        std::cout << "[CODEBOOK INFO] cluster - " << k << std::endl;
        std::cout << "[CODEBOOK INFO] codebook.cluster_locs = " << global.codebook.cluster_locs.row(k) << std::endl;
        std::cout << "[CODEBOOK INFO] codebook.cluster_std = " << global.codebook.cluster_std.row(k)<< std::endl;    // std::cout << "imshowing starting 2" << std::endl;
        std::cout << "[CODEBOOK INFO] codebook.cluster_Ex2 = " << global.codebook.cluster_Ex2.row(k)<< std::endl;    // std::cout << "imshowing starting 2" << std::endl;
        std::cout << "[CODEBOOK INFO] codebook.num_samples = " << global.codebook.n_samples_per_cluster[k]<< std::endl;    // std::cout << "imshowing starting 2" << std::endl;
        
    }
}



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
    cv::reduce(_W, _D, 0, CV_REDUCE_SUM, CV_32F);

    // normalized graph laplacian. See the literature   
    _L =  cv::Mat::eye(_W.size(), CV_32F) - cv::Mat::diag(_D).inv() * _W; 

    
    cv::eigen(_L, _L_eigenval);
    cv::sort(_L_eigenval, _L_eigenval, CV_SORT_EVERY_COLUMN+CV_SORT_ASCENDING);
    cv::sortIdx( _L_eigenval(cv::Range(1,_W.rows),cv::Range::all())-_L_eigenval(cv::Range(0,_W.rows-1),cv::Range::all()), 
        _aux_mat, CV_SORT_EVERY_COLUMN+CV_SORT_DESCENDING);
    n_clusters = _aux_mat.at<int>(0)+1;    

    return n_clusters;
}



void jointsTemporalSmoothingKeyframeSelection() 
// Input: output_vec_KFrameBuffer
// Output: output_vec_KFramePool
// Work: (1) temporal interpolation and (2) keyframe selection
// Notice: this function was expired and is only preserved the functionality for data transfer.!!!
{

    for (int i = 0; i < global.output_vec_KFrameBuffer.size(); i++) {
        global.output_vec_KFramePool.push_back(global.output_vec_KFrameBuffer[i]);
    }
    
    global.output_vec_KFrameBuffer.clear();


}




void initCodebook ( ) {
    
    if (global.VERBOSE)
        std::cout << "[CODEBOOK INFO INIT] compute from data..." << std::endl;
    
    // 1. data standardization
    cv::Mat samples_o = cv::Mat::zeros(global.output_vec_KFramePool.size(), global.num_dims, CV_32F); // the detection score is included
    cv::Mat samples;
    for (int i = 0; i < global.output_vec_KFramePool.size(); i++) 
       global.output_vec_KFramePool[i].copyTo(samples_o.row(i));
       
    
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
    
        
    cv::Mat labels, cluster_locs;
    
    if (global.codebook.init_k_for_kmeans > 1){
        cv::kmeans(samples, global.codebook.init_k_for_kmeans, labels, cv::TermCriteria(cv::TermCriteria::COUNT+cv::TermCriteria::EPS, 1000,1e-6),
               1, cv::KMEANS_PP_CENTERS, cluster_locs);
    }
    else{
        labels = cv::Mat::zeros(samples.rows, 1, CV_32F);
        cv::reduce(samples, cluster_locs, 0, CV_REDUCE_AVG);
    }
    labels.convertTo(labels, CV_32F);

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
    

    if (global.VERBOSE){
        std::cout <<"--compute cluster radius finish" << std::endl;
    }

    if (global.VERBOSE){
        std::cout << "[CODEBOOK INFO INIT] codebook.n_clusters = " << global.codebook.n_clusters << std::endl;
        std::cout << "[CODEBOOK INFO INIT] codebook.cluster_labels = " << global.codebook.cluster_labels.t() << std::endl;
        std::cout << "[CODEBOOK INFO INIT] codebook.cluster_locs = " << global.codebook.cluster_locs << std::endl;
        std::cout << "[CODEBOOK INFO INIT] codebook.cluster_std = " << global.codebook.cluster_std<< std::endl;
        cv::Mat n_samples_per_cluster = cv::Mat(global.codebook.n_samples_per_cluster);
        std::cout << "[CODEBOOK INFO INIT] codebook.n_samples_per_cluster = " << n_samples_per_cluster.t()<< std::endl;
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
        dist.push_back(_dist);
    }
}




void updateCodebook() {
    for(int i = 0; i < global.output_vec_KFramePool.size(); i++) {
        cv::Mat _dist, _joints, _joints_t, _joints_repmat;
        double _dist_min;
        int _dist_min_idx; 
        
        std::vector<float> _cluster_radius;
        _joints = global.output_vec_KFramePool[i];        
        _joints_t = (_joints - global.codebook.init_standardization_mean) 
                        / global.codebook.init_standardization_std;

        calDistSampleToCluster(_joints_t, global.codebook.cluster_locs, _dist);
        cv::minMaxIdx(_dist, &_dist_min, NULL, &_dist_min_idx);
        calClusterRadiusfromStd(_cluster_radius);
        // float _radius_max = *std::max_element(_cluster_radius.begin(), _cluster_radius.end() );
        

        if (global.VERBOSE){
            std::cout << "[DEBUG] _dist_min = " << _dist_min << "  [DEBUG] _radius_min_dist = " << _cluster_radius[_dist_min_idx] << std::endl;
        }

        if ( _dist_min > _cluster_radius[_dist_min_idx]) { 

            if (global.VERBOSE)
                std::cout<< "[CODEBOOK INFO] updateCodebook(): creating new cluster";

            global.codebook.cluster_labels.push_back(global.codebook.n_clusters);
            global.codebook.sample_labels.push_back((float)global.codebook.n_clusters);
            global.codebook.n_clusters += 1;
            global.codebook.cluster_locs.push_back(_joints_t);

            
            cv::Mat zz = cv::Mat::zeros(1,_joints_t.cols, CV_32F);
            cv::Mat _ex2 = _joints_t.mul(_joints_t);
            global.codebook.cluster_std.push_back(zz);
            global.codebook.cluster_Ex2.push_back(_ex2);
            global.codebook.n_samples_per_cluster.push_back(1);


        }
        else {  
            
            
            if (global.VERBOSE)
                std::cout<< "[CODEBOOK INFO] updateCodebook(): updating current cluster " << _dist_min_idx << std::endl;
            
            global.codebook.sample_labels.push_back((float)_dist_min_idx);
            int _n_samplei = global.codebook.n_samples_per_cluster[_dist_min_idx];
            cv::Mat _mean = global.codebook.cluster_locs.row(_dist_min_idx).clone();
            cv::Mat _std = global.codebook.cluster_std.row(_dist_min_idx).clone();
            cv::Mat _Ex2 = global.codebook.cluster_Ex2.row(_dist_min_idx).clone();
            

            // updating mean
            global.codebook.cluster_locs.row(_dist_min_idx) = _mean + (_joints_t-_mean ) / (float) (1+_n_samplei);

            // updating Ex2
            global.codebook.cluster_Ex2.row(_dist_min_idx) = (_Ex2*_n_samplei + _joints_t.mul(_joints_t))/(float)(_n_samplei+1);

            // updating std
            _std = global.codebook.cluster_Ex2.row(_dist_min_idx) 
                  - global.codebook.cluster_locs.row(_dist_min_idx).mul(global.codebook.cluster_locs.row(_dist_min_idx));

            cv::sqrt(_std, global.codebook.cluster_std.row(_dist_min_idx));        

            // update sample numbers
            global.codebook.n_samples_per_cluster[_dist_min_idx] = _n_samplei + 1;
            
        }

    }
}








// in this function, keyframe candidate selection and incremental clustering run successively
void runIncrementalClustering() {

    jointsTemporalSmoothingKeyframeSelection();

    if (global.VERBOSE){
        std::cout << "[CODEBOOK INFO] global.KEYFRAME_POOL_FILLTIMES=" << global.KEYFRAME_POOL_FILLTIMES << std::endl;
        std::cout << "[CODEBOOK INFO] global.KEYFRAME_POOL_SIZE=" << global.output_vec_KFramePool.size() << std::endl;
    }
    
    if (global.KEYFRAME_POOL_FILLTIMES ==0) { 
        if (global.VERBOSE)
            std::cout << "[CODEBOOK INFO] runIncrementalClglobal.n_featurestering(): Codebook Initialization--------------------BEGIN" << std::endl;
        initCodebook( );

        if (global.VERBOSE)
            std::cout << "[CODEBOOK INFO] runIncrementalClustering(): Codebook Initialization--------------------END" << std::endl;

    }
    else {
        if (global.VERBOSE)
            std::cout << "[CODEBOOK INFO] runIncrementalClustering(): Codebook Updating--------------------BEGIN" << std::endl;

            updateCodebook( );

        if (global.VERBOSE)
            std::cout << "[CODEBOOK INFO] runIncrementalClustering(): Codebook Updating--------------------END" << std::endl;
    }
    
    global.KEYFRAME_POOL_FILLTIMES+=1;
    global.output_vec_KFramePool.clear();

}



// Input: all the pose patterns. Rows - samples, cols - features
// Output: codebook, segment boudnaries.
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray *prhs[])
{

    if (nrhs != 4)
        mexErrMsgIdAndTxt("Three arguments should be parsed", "usage: fun(samples, ShortMemorySize, errorEstimate, VERBOSE)");


    // parse arguments
    global.KEYFRAME_BUFFERSIZE = (int)*mxGetPr(prhs[1]);
    int offset = global.KEYFRAME_BUFFERSIZE;
    float error_estimate = (float)*mxGetPr(prhs[2]);
    
    global.VERBOSE = (int)*mxGetPr(prhs[3]);
    
    cv::Mat all_samples;
    ocvMxArrayToMat_double(prhs[0], all_samples);
    all_samples.convertTo(all_samples, CV_32F);
    
    global.num_samples = all_samples.rows;
    global.num_dims = all_samples.cols;
    global.KEYFRAME_POOL_FILLTIMES = 0;
    global.delta = error_estimate*(global.num_dims);
    global.output_vec_KFrameBuffer.clear();
    global.output_vec_KFramePool.clear();
    global.codebook.cluster_labels.release();
    global.codebook.sample_labels.release();
    global.codebook.cluster_Ex2.release();

    // step2: run the incremental learning mechanism
    for(int i = 0; i < global.num_samples; i+=offset){
        
        int upper_idx_bound = std::min(i+global.KEYFRAME_BUFFERSIZE, global.num_samples);
        cv::Mat samples_for_buffer = all_samples.rowRange(i, upper_idx_bound);
        for(int j = 0; j < samples_for_buffer.rows; j++)
            global.output_vec_KFrameBuffer.push_back(samples_for_buffer.row(j));

        runIncrementalClustering();
    }
    if(global.VERBOSE)
        showCodebookInfo();

    cv::Mat labels;
    cv::Mat cluster_locs;
    global.codebook.sample_labels.convertTo(labels, CV_64F);
    global.codebook.cluster_locs.convertTo(cluster_locs, CV_64F);
    if (global.VERBOSE)
        std::cout << "before write to mxArray" << std::endl;

    plhs[0] = ocvMxArrayFromMat_double(labels);
    plhs[1] = ocvMxArrayFromMat_double(cluster_locs);


}

