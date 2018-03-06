function idx = calLocalFeatureAggregationAndClustering(X,Xl,C,n_clusters,varargin)


% Xl(Xl~=0) = 1;
%%% find the peaks in the label space to determine the time window
% labels = imgaussfilt(Xl,12.5); %% for CMUMAD
labels = imgaussfilt(Xl,1.5); %% for TUMKitchen
% peak_width_weight = 1; %% for CMUMAD
peak_width_weight = 0.5; %% for TUMKitchen

[pks_label,locs_label, pks_width] = findpeaks(labels);
% pks_width = 2*pks_width;
pks_width = round(pks_width);
%%% use the time window to aggregate features
features = zeros(length(pks_label), size(C,1));
for pp = 1:length(pks_label)
    lb = max(1,locs_label(pp)-round(peak_width_weight*pks_width(pp)));
    ub = min(size(X,1), locs_label(pp)+round(peak_width_weight*pks_width(pp)));
    XX = X(lb : ub, :);
    dist = pdist2(XX, C);
    ff = exp(-0.1*dist)./ repmat(sum(exp(-0.1*dist),2), 1, size(C,1));
    features(pp,:) = sum(ff,1);
    features(pp,:) = features(pp,:)/norm(features(pp,:));
end


method = 'Kmeans';
if nargin == 4
    method = 'Kmeans';
    nc = n_clusters-1;
else
    agg_type = varargin{1};
end

if strcmp(agg_type, 'null_action')
    nc = n_clusters-1;
elseif strcmp(agg_type, 'normal_action')
    nc = n_clusters-1;
end


if size(features,1) <= nc
    iidx = 1:size(features,1);
elseif strcmp(method, 'Kmeans')
    iidx = kmeans(features, nc);
    
elseif strcmp(method,'ours')
    time_window = 50;
    sigma = 0.1;
    
    [iidx, ~] = dynamic_clustering(features, time_window,sigma,0);

end


idx = zeros(size(Xl));
for pp = 1:length(pks_label)
    lb = max(1,locs_label(pp)-pks_width(pp));
    ub = min(size(X,1), locs_label(pp)+pks_width(pp));
    idx(lb:ub)=iidx(pp);
end

    

