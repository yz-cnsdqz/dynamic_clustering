function idx = calLocalFeatureAggregationAndClustering_TUMKitchen(X,Xl,C,n_clusters,varargin)


%%% sliding window method
W = 30;S = 1; 
n_frames = size(X,1);


%%% use the time window to aggregate features
features = zeros(n_frames,size(C,1));
for ii = 1:S:n_frames
    lb = max(1,ii-W/2);
    ub = min(n_frames,ii+W/2);
    XX = X(lb : ub, :);
    dist = pdist2(XX, C);
    ff = exp(-0.1*dist)./ repmat(sum(exp(-0.1*dist),2), 1, size(C,1));
    features(ii,:) = sum(ff,1);
    features(ii,:) = features(ii,:)/norm(features(ii,:));
end


idx = kmeans(features, n_clusters);

end

    

