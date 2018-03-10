function idx = calLocalFeatureAggregationWithEncodedFeatures_TUMKitchen(X,Xl,n_clusters,varargin)


%%% sliding window method
W = 30;S = 1; 
n_frames = size(X,1);


%%% use the time window to aggregate features
features = zeros(n_frames,size(X,2));
for ii = 1:S:n_frames
    lb = max(1,ii-W/2);
    ub = min(n_frames,ii+W/2);
    XX = X(lb : ub, :);
    features(ii,:) = sum(XX,1);
    features(ii,:) = features(ii,:)/norm(features(ii,:));
end


idx = kmeans(features, n_clusters);

end

    

