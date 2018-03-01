clear all;
close all;
clc;
addpath(genpath('../../eval_package'));
addpath(genpath('../../aca'));
addpath(genpath('../../TSC'));
addpath(genpath('../../mexIncrementalClustering'));
dataset_path = '/mnt/hdd/Dataset_TUMKitchen2';
video_list = importdata([dataset_path '/video_list.txt']);


%%% scenario configuration - before each running, check here!"
feature_list = {'jointLocs','relativeAngle','quaternion'};
method_list = {'kmeans','spectralClustering','TSC','ACA','ours'};
bodypart_list = {'rightArm','leftArm','torso'};

method = method_list{end};
bodypart = bodypart_list{3};
feature = feature_list{1};

is_show = 1; 
%%% scenario configuration - end"

Precision = [];
Recall = [];
CptTime = [];
% for vv = 1:length(video_list)
for vv = 10:10    
    %------------------ fearture extraction and annotation----------------%
    skeleton = importdata([dataset_path '/' video_list{vv} '/poses.csv']);    
    [xt, xl, xr] = calPatternFromSkeleton(skeleton, feature);
    anno = importdata([dataset_path '/' video_list{vv} '/labels.csv']);
    [yt_left_arm, yt_right_arm, yt_torso] = readTUMKitchenAnnotation(anno);
    if strcmp(bodypart,'leftArm')
        n_clusters = 9;
        yt = yt_left_arm;
        pattern = xl;
    elseif strcmp(bodypart,'rightArm')
        n_clusters = 9;
        yt = yt_right_arm;
        pattern = xr;
    elseif strcmp(bodypart,'torso')
        n_clusters = 2;
        yt = yt_torso;
        pattern = xt;
    end
    %------------------ run methods and make evaluations------------------%
    
    startTime = tic;
    if strcmp(method, 'kmeans')
%         [idx,C] = kmeans(pattern,n_clusters);
        [idx,C] = kmeans(pattern,50);
        
%         idx = calLocalFeatureAggregationAndClustering_TUMKitchen(pattern,idx,C, n_clusters,'normal_action'); 

    elseif strcmp(method, 'spectralClustering')
        idx = spectralClustering(pattern,500,n_clusters);
        idx = idx-1; %%% 0-based label

    elseif strcmp(method, 'TSC')
        %%%---Normalize the data---%%%
%             X = normalize(pattern);

        %%%---Parameter settings---%%%
        paras = [];
        paras.lambda1 = 0.01;
        paras.lambda2 = 15;
        paras.n_d = 80;
        paras.ksize = 7;
        paras.tol = 1e-4;
        paras.maxIter = 12;
        paras.stepsize = 0.1;

        %%%---Learn representations Z---%%%
%             disp('--first pca to 100d; otherwise computation is prohibitively expensive.');
%             [comp,XX,~] = pca(pattern, 'NumComponents',100);
        [D, Z, err] = TSC_ADMM(pattern',paras);
        disp('clustering via graph cut..');
%             nbCluster = length(unique(label));
        vecNorm = sqrt(sum(Z.^2));
        W2 = (Z'*Z) ./ (vecNorm'*vecNorm + 1e-6);
        [oscclusters,~,~] = ncutW(W2,n_clusters);
        idx = denseSeg(oscclusters, 1);
        idx = idx;

%         uid = idx(1);
%         idx(idx==uid) = 10e6;
%         idx(idx==1) = uid;
%         idx(idx==10e6) = 1;
        idx = idx-1; %%% 0-based label


    elseif strcmp(method, 'ACA')
        idx = calACAOrHACA(pattern,n_clusters, 'ACA');
        idx(end) = []; %%% remove redudant frame

%         uid = idx(1);
%         idx(idx==uid) = 10e6;
%         idx(idx==1) = uid;
%         idx(idx==10e6) = 1;
%         idx = idx-1; %%% 0-based label
%         save('ACA_idx_TUMKitchen.mat','idx');

    elseif strcmp(method, 'HACA')
        idx = calACAOrHACA(pattern,36, 'HACA');
        idx(end) = []; %%% remove redudant frame

    elseif strcmp(method,'ours')
        time_window = 50;
        sigma = 25;

%         disp('--online learn the clusters and labels..');
        [idx1, C] = incrementalClustering(double(pattern), time_window,sigma,0,0,1.0);
%         idx2 = zeros(size(idx1));
%         for kk = 2:length(idx1)
%             if idx1(kk)~=idx1(kk-1)
%                 idx2(kk) = idx1(kk)+randi(length(unique(idx1)))-1;
%             else
%                 idx2(kk) = idx2(kk-1);
%             end
%         end
%         uid = mode(idx2);
%         idx2(idx2==uid) = 10e6;
%         idx2(idx2==0) = uid;
%         idx2(idx2==10e6) = 0;
%         disp('--postprocessing, merge clusters');
%         idx = calLocalFeatureAggregationAndClustering_TUMKitchen(pattern,idx1,C, n_clusters,'normal_action'); 
        idx = idx1;
    end
    CptTime = [CptTime toc(startTime)];
%     Result= funEvalSegmentation_TUMKitchen(yt, idx, 7, is_show, method,bodypart,feature,vv);
    Result= funEvalSegmentation_TUMKitchen(yt, idx, 7, is_show);
    Precision = [Precision Result.Prec];
    Recall = [Recall Result.Rec];
end
disp('====================================================')
disp('Evaluation Results:')
fprintf('Method: %s\n',method);
fprintf('Bodypart: %s\n',bodypart);
fprintf('Feature: %s\n',feature);
fprintf('ave_precision: %f\n',mean(Precision));
fprintf('ave_recall: %f\n',mean(Recall));
fprintf('ave_runtime: %f\n',mean(CptTime));



    
    
