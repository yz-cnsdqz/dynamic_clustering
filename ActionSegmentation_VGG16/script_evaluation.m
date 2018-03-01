clear all;
close all;
clc;
addpath(genpath('../../eval_package'));
addpath(genpath('../../TSC'));
addpath(genpath('../../mexIncrementalClustering'));

%% evaluate script
dataset = '/home/yzhang/Videos/Dataset_CMUMAD';
subject = 1:20;
feature = 'VGG';
% method_set = {'TSC','spectralClustering', 'ACA'};
method_set = {'kmeans'};


for mm = 1:length(method_set)
    method = method_set{mm};

    fprintf('======== use %s ===============\n', method);


    Pre = [];
    Rec = [];
    Acc = [];
    CptTime = [];
    for ss = subject
        pre = 0;
        rec = 0;
        acc = 0;    
        for qq = 1:2

            feature_file = sprintf([dataset, '/',feature,'/CMUMAD_sub%02d_seq%02d.mat'], ss,qq);
            data = load(feature_file);
            gt_file = sprintf([dataset,'/sub%02d/seq%02d_label.mat'], ss,qq);
            load(gt_file); %  subject 1, sequence 
            gtlabE=  extractGTFormat(label); % true event-based labels (obtain by loading a true label file)
            n_frames = label(end,3)-1;
            pattern1 = double(data.pattern');
    %         pattern = pattern1(:,1:101)+1.0*pattern1(:,102:end);
            pattern = [pattern1; pattern1(end,:)];

    %         pattern = prdInterpolation(pattern, 50, n_frames); % stride and time window is fixed.

            n_clusters = 36;
            startTime = tic;
            if strcmp(method, 'kmeans')
                [idx,C] = kmeans(pattern,n_clusters);
                uid = idx(1);
                idx(idx==uid) = 10e6;
                idx(idx==1) = uid;
                idx(idx==10e6) = 1;
                idx = idx-1; %%% 0-based label
                idx = calLocalFeatureAggregationAndClustering(pattern,idx,C, 36); 


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
                disp('--first pca to 100d; otherwise computation is prohibitively expensive.');
                [comp,XX,~] = pca(pattern, 'NumComponents',100);
                [D, Z, err] = TSC_ADMM(XX',paras);
                disp('clustering via graph cut..');
    %             nbCluster = length(unique(label));
                vecNorm = sum(Z.^2);
                W2 = (Z'*Z) ./ (vecNorm'*vecNorm + 1e-6);
                [oscclusters,~,~] = ncutW(W2,36);
                idx = denseSeg(oscclusters, 1);
                idx = idx';







            elseif strcmp(method, 'ACA')
                addpath(genpath('../../aca'));

                idx = calACAOrHACA(pattern,36, 'ACA');
                idx(end) = []; %%% remove redudant frame
                rmpath(genpath('../../aca'));



            elseif strcmp(method, 'HACA')
                idx = calACAOrHACA(pattern,36, 'HACA');
                idx(end) = []; %%% remove redudant frame




            elseif strcmp(method, 'KTC-S')
                disp('--first pca to 100d; otherwise computation is prohibitively expensive.');
                [comp,XX,~] = pca(pattern, 'NumComponents',100);
                sigma_s = 0.001;
                sigma_t = 0.001;
                time_window = 50;
                idx = calKernelizedTemporalSegment(XX,time_window,sigma_s,sigma_t);


            elseif strcmp(method,'ours')
                time_window = 50;
                sigma = 1;
                is_temporal_reg = 0;

                disp('--online learn the clusters and labels..');
                [idx, C] = incrementalClustering(double(pattern), time_window,sigma,is_temporal_reg,0,2.5);
                %disp('--postprocessing, merge clusters');
                idx = calLocalFeatureAggregationAndClustering(pattern,idx,C, 36); 

            end
            CptTime = [CptTime toc(startTime)];
            res = funEvalDetection(gtlabE, double(idx'), 0.5, false); % tol range = 0.5
            %pre = pre+res.Prec/2;
            %rec = rec+res.Rec/2;
            %acc = acc+res.Acc/2;
            Pre = [Pre res.Prec];
            Rec = [Rec res.Rec];
            Acc = [Acc res.Acc];

        end
        fprintf('--subject %02d: precision=%f, recall=%f, accuracy=%f\n',ss,pre,rec,acc);
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
    
    result_filename = sprintf('Result_VGG_%s.mat', method);
    save(result_filename, 'Result');
end



