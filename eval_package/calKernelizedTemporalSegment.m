function labels = calKernelizedTemporalSegment(X, T, sigma_s, sigma_t)
%%% this script implement KTC-S methed introduced in Dian Gong, et al. 2012
local_set = 5; %%% when compute local tagent similarity, use this local set.
[n_samples, n_dims] = size(X);
labels = zeros(n_samples,1);

%%% main loop
p = 1;
ll = 0;
while(p < n_samples)
    up = min(p+T-1, n_samples);
    Xt = X(p:up, :);
%     disp('-- compute kernel');
    Kst = calSpatialTemporalKernelMatrix(Xt, sigma_s, sigma_t, local_set);
%     Kst = real(Kst);
    %%% optimize the energy function via searching
%     disp('-- optimize..');
    cut = 2;
    loss_min = 10^10;
    e0 = zeros(T,1);
    up2 = min(T, size(Kst,2));
    for ii = 2:up2-1
        k11 = Kst(1:ii, 1:ii);
        k22 = Kst(ii+1:up2, ii+1:up2);
        k12 = Kst(ii+1:up2, 1:ii);
        loss = sum(k11(:))/prod(size(k11)) + sum(k22(:))/prod(size(k22)) ...
            - 2*sum(k12(:))/prod(size(k12));
        if loss < loss_min
            cut = ii;
            loss_min = loss;
        end
    end
    cut = p+cut;
    %%% assign labels depending on the cut
    ll = ll+1;
    labels(p:cut) = ll;
    ll = ll+1;
    labels(cut+1:up) = ll;
    
    p = cut+1;
%     fprintf('-- cut = %d\n', cut);
end

labels(labels==0) = ll;

end



function Kst = calSpatialTemporalKernelMatrix(X, sigma_s, sigma_t, local_set)

[n_samples, n_dims] = size(X);

%%% spatial kernel matrix;
Ks = exp(-sigma_s * squareform(pdist(X)));

%%% temporal kernel matrix;
% V = calLocalTangentNormal(X, local_set)-1e-6;
% Kt= exp(-sigma_t* acos(V*V'));

%%% combine them
% Kst = Ks.*Kt;
Kst = Ks;
end


function V = calLocalTangentNormal(X, local_set)

n_samples = size(X,1);
V = [];
%%% mirror the boundary
bd = (local_set-1)/2;
XX = zeros(size(X,1)+2*bd, size(X,2));
XX(bd+1 : bd+n_samples,:) = X;
for ii = 1:bd
    XX(ii,:) = X(1,:);
    XX(end+1-ii,:) = X(end,:);
end

%%% calculate tangent
for ii = 1:n_samples
    [comp] = pca(XX(ii : ii+2*bd, :)'*XX(ii : ii+2*bd, :));
    V = [V; comp(:,end)'];
end

end
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        



