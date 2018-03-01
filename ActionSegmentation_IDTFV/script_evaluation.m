clear all;
close all;
clc;
addpath(genpath('../../eval_package'));
addpath(genpath('../../TSC'));
addpath(genpath('../../mexIncrementalClustering'));

%% evaluate script
dataset = '/home/yzhang/Videos/Dataset_CMUMAD';
%  subject = 1:20;
subject = 5;
feature = 'FV2';
method_set = {'ACA'};


for mm = 1:length(method_set)
    method = method_set{mm};
    fprintf('============= evaluation method: %s==================\n',method);

    Pre = [];
    Rec = [];
    Acc = [];
    CptTime = [];

    for ss = subject
        pre = 0;
        rec = 0;
        acc = 0;
        feature_file = sprintf([dataset, '/',feature,'/CMUMAD_EvaluationResults_sub%d.mat'], ss);
        load(feature_file);
        
        for qq = 1:2
            if qq ==1
                qqs = 2;
            else
                qqs = 1;
            end
            gt_file = sprintf([dataset,'/sub%02d/seq%02d_label.mat'], ss,qq);
            load(gt_file); %  subject 1, sequence 
            gtlabE=  extractGTFormat(label); % true event-based labels (obtain by loading a true label file)
            n_frames = label(end,3);

            pattern=eval_res{qqs}.stip_T_encoded{1}.feature; % test frame-based labels produced by SVM+DP (baseline)
            
            pattern = prdInterpolation(pattern, 50, n_frames); % stride and time window is fixed.
    %         pattern = [pattern; pattern(end,:)];
            n_clusters = 36;
            startTime = tic;
            if strcmp(method, 'kmeans')
                [idx,C] = kmeans(pattern,n_clusters);
                uid = idx(1);
                idx(idx==uid) = 10e6;
                idx(idx==1) = uid;
                idx(idx==10e6) = 1;
                idx = idx-1; %%% 0-based label
    %             idx = calLocalFeatureAggregationAndClustering(pattern,idx,C, 36); 




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

                uid = idx(1);
                idx(idx==uid) = 10e6;
                idx(idx==1) = uid;
                idx(idx==10e6) = 1;
                idx = idx-1; %%% 0-based label





            elseif strcmp(method, 'ACA')
                [comp,XX,~] = pca(pattern, 'NumComponents',100);
                addpath(genpath('../../aca'));

                idx = calACAOrHACA(XX,36, 'ACA');
                idx(end) = []; %%% remove redudant frame
                rmpath(genpath('../../aca'));
                uid = idx(1);
                idx(idx==uid) = 10e6;
                idx(idx==1) = uid;
                idx(idx==10e6) = 1;
                idx = idx-1; %%% 0-based label
                rmpath(genpath('../../aca'));
 



            elseif strcmp(method, 'HACA')
                [comp,XX,~] = pca(pattern, 'NumComponents',100);
                addpath(genpath('../../aca'));
                idx = calACAOrHACA(pattern,36, 'HACA');
                idx(end) = []; %%% remove redudant frame
                rmpath(genpath('../../aca'));




            elseif strcmp(method, 'KTC-S')
                disp('--first pca to 100d; otherwise computation is prohibitively expensive.');
                [comp,XX,~] = pca(pattern, 'NumComponents',100);
                sigma_s = 0.001;
                sigma_t = 0.001;
                time_window = 50;
                idx = calKernelizedTemporalSegment(XX,time_window,sigma_s,sigma_t);


            elseif strcmp(method,'ours')
                time_window = 50;
                sigma = 0.0075;

    %             disp('--online learn the clusters and labels..');
                [idx, C] = incrementalClustering(double(pattern), time_window,sigma,0, 0 ,1.0);

    %             disp('--postprocessing, merge clusters');
    %             idx = calLocalFeatureAggregationAndClustering(pattern,idx,C, 36); 

            end
            CptTime = [CptTime toc(startTime)];
            res = funEvalDetection(gtlabE, double(idx'), 0.5, true); % tol range = 15
    %         pre = pre+res.Prec/2;
    %         rec = rec+res.Rec/2;
    %         acc = acc+res.Acc/2;
            Pre = [Pre res.Prec];
            Rec = [Rec res.Rec];
            Acc = [Acc res.Acc];
        end

        fprintf('--subject %02d: finish..\n',ss);
    end

    fprintf('- avg precision = %f\n', mean(Pre));
    fprintf('- avg recall = %f\n', mean(Rec));
    fprintf('- avg accuracy = %f\n', mean(Acc));
    fprintf('- avg runtime = %f\n', mean(CptTime));
    
    
    Result.Precision = Pre;
    Result.Recall = Rec;
    Result.Accuracy = Acc;
    Result.RunTime = CptTime;

    Result.meanPrecision = mean(Pre);
    Result.meanRecall = mean(Rec);
    Result.meanAccuracy = mean(Acc);
    Result.meanRunTime = mean(CptTime);
    result_filename = sprintf('Result_FV_%s.mat',method);
    save(result_filename, 'Result');
    
    
end


