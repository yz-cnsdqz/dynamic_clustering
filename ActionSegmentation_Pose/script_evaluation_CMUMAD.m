clear all;
close all;
clc;
addpath(genpath('../../eval_package'));
addpath(genpath('../../aca'));
addpath(genpath('../../TSC'));
addpath(genpath('../../mexIncrementalClustering'));
%% evaluate script
dataset = '/home/yzhang/Videos/Dataset_CMUMAD';
% subject = 1:20;
subject = 1:20;
feature = 'BodyPose';
%method_set = {'spectralClustering', 'TSC','ACA'};
%pose_feature_set = {'jointLocs'};
method_set = {'ours'};
pose_feature_set = {'jointLocs', 'relativeAngle', 'quaternions'};
time_weights = 0.3*1e-2;

for mm = 1:length(method_set)
    for pp = 1:length(pose_feature_set)
        
        method = method_set{mm};
        pose_feature = pose_feature_set{pp};
        fprintf('================ perform %s ===== %s =========\n', method, pose_feature);


        Pre = [];
        Rec = [];
        Acc = [];
        CptTime = [];
        for ss = subject
            pre = 0;
            rec = 0;
            acc = 0;

            for qq = 1:2
                feature_file = sprintf([dataset, '/',feature,'/PoseFeature_sub%02d_seq%02d.mat'], ss,qq);
                data = load(feature_file);
                gt_file = sprintf([dataset,'/sub%02d/seq%02d_label.mat'], ss,qq);
                load(gt_file); %  subject 1, sequence   
                gtlabE =  extractGTFormat(label); % true event-based labels (obtain by loading a true label file)
                n_frames = label(end,3)-1;

                if strcmp(pose_feature, 'jointLocs')
                    pattern = data.pattern.jointLocs;
                elseif strcmp(pose_feature,'relativeAngle')
                    pattern = data.pattern.relativeAngle;
                elseif strcmp(pose_feature, 'quaternions')
                    pattern = data.pattern.quaternions;
                end

                n_clusters = 36;
                startTime = tic;
                if strcmp(method, 'kmeans')
                    [idx,C] = kmeans(pattern,36);
                    uid = idx(1);
                    idx(idx==uid) = 10e6;
                    idx(idx==1) = uid;
                    idx(idx==10e6) = 1;
                    idx = idx-1; %%% 0-based label
                    idx = calLocalFeatureAggregationAndClustering(pattern,idx,C, 36); 

                    
                elseif strcmp(method, 'regularSampling')
                        %%% regular sampling
                    prd_boundary = 1:30:size(pattern,1);
                    idx = zeros(size(pattern,1),1);
                    idx(prd_boundary) = 1;
                    
                elseif strcmp(method, 'spectralClustering')
                    idx = spectralClustering(pattern,1,n_clusters);
                    uid = idx(1);
                    idx(idx==uid) = 10e6;
                    idx(idx==1) = uid;
                    idx(idx==10e6) = 1;
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
                    vecNorm = sum(Z.^2);
                    W2 = (Z'*Z) ./ (vecNorm'*vecNorm + 1e-6);
                    [oscclusters,~,~] = ncutW(W2,36);
                    idx = denseSeg(oscclusters, 1);
                    idx = idx';

                    uid = idx(1);
                    idx(idx==uid) = 10e6;
                    idx(idx==1) = uid;
                    idx(idx==10e6) = 1;
                    idx = idx-1; %%% 0-based label



                elseif strcmp(method, 'KTC-S')
        %             disp('--first pca to 100d; otherwise computation is prohibitively expensive.');
        %             [comp,XX,~] = pca(pattern, 'NumComponents',100);
                    sigma_s = 1;
                    sigma_t = 0.001;
                    time_window = 50;
                    idx = calKernelizedTemporalSegment(pattern,time_window,sigma_s,sigma_t);


                elseif strcmp(method, 'ACA')
                    idx = calACAOrHACA(pattern,36, 'ACA');
                    idx(end) = []; %%% remove redudant frame
                    
                    uid = idx(1);
                    idx(idx==uid) = 10e6;
                    idx(idx==1) = uid;
                    idx(idx==10e6) = 1;
                    idx = idx-1; %%% 0-based label
                    save('ACA_idx.mat','idx');

                elseif strcmp(method, 'HACA')
                    idx = calACAOrHACA(pattern,36, 'HACA');
                    idx(end) = []; %%% remove redudant frame

                elseif strcmp(method,'ours')
                    time_window = 10;
                    sigma = 0.01;

%                     disp('--online learn the clusters and labels..');
                    [idx, C] = incrementalClustering(double(pattern), time_window,sigma,0);

                    %%% uncomment the following for online processing
%                     disp('--postprocessing, merge clusters');
%                     idx = calLocalFeatureAggregationAndClustering(pattern,idx,C, 36); 

                end
                CptTime = [CptTime toc(startTime)];
%                 res = funEvalDetection(gtlabE, double(idx'), 0.5, true); % tol range = 15
                yt = zeros(size(idx));
                yt(cumsum(gtlabE(:,2)))=1;
                yt(length(idx)+1:end) = [];
                res= funEvalSegmentation_TUMKitchen(yt, idx, 7, 0);

                %pre = pre+res.Prec/2;
                %rec = rec+res.Rec/2;
                %acc = acc+res.Acc/2;
                Pre = [Pre res.Prec];
                Rec = [Rec res.Rec];
%     	        Acc = [Acc res.Acc];

            end
            
            %fprintf('--subject %02d: precision=%f, recall=%f, accuracy=%f\n',ss,pre,rec,acc);
        end

        fprintf('- avg_pre = %f\n',mean(Pre));
        fprintf('- avg_rec = %f\n',mean(Rec));
        fprintf('- avg_acc = %f\n',mean(Acc));
        fprintf('- avg_runtime = %f seconds\n',mean(CptTime));
        Result.Precision = Pre;
        Result.Recall = Rec;
        Result.Accuracy = Acc;
        Result.RunTime = CptTime;

        Result.meanPrecision = mean(Pre);
        Result.meanRecall = mean(Rec);
        Result.meanAccuracy = mean(Acc);
        Result.meanRunTime = mean(CptTime);
        result_filename = sprintf('Result_%s_%s.mat',pose_feature,method);
        save(result_filename, 'Result');
    end
end
