import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from sklearn.cluster import KMeans



def compute_distance_to_clusters(x, centers, dist_type):
    ''' Identify closest center & compute distance to it '''
    centers_array = np.concatenate(centers,axis=0)
    dist_array = cdist(x, centers_array, metric=dist_type) #compute pairwise distance between x & centers
    min_dist = np.amin(dist_array)
    min_idx = np.argmin(dist_array)
    return min_dist, min_idx


def graph_spectral_init(inputs, gamma=1.0):
    
    centers = [] #cluster centers
    radius = [] #cluster radius
    idx = []    # cluster index
    M = []      # helper variable
    NSamples = []   # cluster size
    
    # first, create a fully connected graph with gaussian similarity
    X = np.concatenate(inputs, axis=0)
    n_samples, n_dims = X.shape
    A = squareform(pdist(X, metric='sqeuclidean'))
    W = np.exp(-gamma*A) - np.identity(n_samples)
    D = np.sum(W, axis=0)


    # second, compute normalized graph laplacian
    L = np.identity(n_samples) - np.matmul(np.linalg.inv(np.diag(D)), W) 
    L_eigenval, L_eigenvec = np.linalg.eig(L)
    sort_idx = np.argsort(L_eigenval)
    L_eigenval = L_eigenval[sort_idx]
    

    # third, select the most significant jump for determining number of clusters
    L_eigenval_diff = L_eigenval[1:]-L_eigenval[:-1]
    n_clusters = np.argmax(L_eigenval_diff)+1
    kmeans = KMeans(n_clusters=n_clusters).fit(X)
    centers= np.split(kmeans.cluster_centers_, n_clusters, axis=0)
    idx = np.split(kmeans.labels_-1, n_samples, axis=0)
    
    for i in range(n_clusters):
        NSamples.append(idx.count(i))
        M.append(X[kmeans.labels_==i+1,:]*X[kmeans.labels_==i+1,:] / NSamples[i])    
        radius.append(M[i]-centers[i]*centers[i])

    return centers, radius, idx, M, NSamples



def dynamic_clustering(inputs, delta, verbose=False, init_with_sc=False, fps=None, graph_weights=1e-5, ClusterStructure=None, dist_type='sqeuclidean'):
## Inputs:
# x:                 input data array, with shape [n_frames, n_dimensions]
# delta:             the hyper-parameter to control the complexity of cluster structure. delta=inf gives only one cluster and delta=0 gives n_frames cluster.
# init_with_sc:      If true, using the first 1-second frames to initialize the cluster structure; else, the first feature vector is regarded as the initial cluster structure
# fps:               video parameter, frames-per-second. Not used if init_with_sc=False
# graph_weights:     the edge weight of the fully-connected graph for initialization. Not used if init_with_sc=False
# cluster_structure: If none, the algorithm runs from scratch. Or the algorithm continues to update the ClusterStructure using the data inputs, and other settings are overwritten by ClusterStructure.



## Outputs: ClusterStructure
# idx:          cluster index of each sample, list
# centers:      cluster centers, list
# radius:       cluster radius, list
# delta:        the delta for updating clusters, number
# NSamples:     cluster size, list
# M:            helper variable, list, E[X^2] for each cluster


    n_frames, n_dims = inputs.shape
    X = np.split(inputs, n_frames, axis=0)
    idx_saver = 0
    ii = 0

    if ClusterStructure is None:        
        if verbose:
            print('[RUNNING LOG] initial cluster structure is empty.')

        centers = [] #cluster centers
        radius = [] #cluster radius
        idx = []    # cluster index
        M = []      # helper variable
        NSamples = []   # cluster size

    else:
        if verbose:
            print('[RUNNING LOG] continue to update the given cluster structure.')

        centers = ClusterStructure['centers']
        radius = ClusterStructure['radius']
        idx = ClusterStructure['idx']
        idx_saver = idx[-1]
        delta = ClusterStructure['delta']
        M = ClusterStructure['M']
        NSamples = ClusterStructure['NSamples']
        init_with_sc = False

    dr = n_dims * delta

    if not init_with_sc:
        if verbose:
            print('[RUNNING LOG] the first sample is regarded as the first cluster.')
        
        for x in X:
            ii += 1
            # if verbose:
            #     print('- iteration={:i}'.format(ii))

            if not centers: # if center is empty
                centers.append(x)
                idx.append(idx_saver)
                radius.append(0.0)
                NSamples.append(1.0)
                M.append(x*x)
            else:

                min_dist, min_idx = compute_distance_to_clusters(x, centers, dist_type)
                # if verbose:
                #     print('-- the min_dist={:f}, min_idx={:d}'.format(min_dist, min_idx))
                #     print('-- the thresh={:f}'.format(dr))

                if min_dist > max(radius[min_idx], dr):
                    centers.append(x)
                    NSamples.append(1.0)
                    M.append(x*x)
                    radius.append(0.0)
                    idx_saver +=1
                    idx.append(idx_saver)
                else:
                    NSamples[min_idx] += 1
                    centers[min_idx] = centers[min_idx]+(x-centers[min_idx])/(NSamples[min_idx])
                    M[min_idx] = M[min_idx] + (x*x - M[min_idx])/(NSamples[min_idx])
                    Sigma = M[min_idx] - centers[min_idx]*centers[min_idx]
                    radius[min_idx] = np.sum(Sigma)
                    idx.append(min_idx)


    else:
        if verbose:
            print('[RUNNING LOG] obtain the initial cluster structure via graph spectra, using the first 1-second features')

        ii = fps
        centers, radius, idx, M, NSamples = graph_spectral_init(X[:fps], graph_weights)    
        idx_saver = idx[-1]
        for x in X[fps:]:

            min_dist, min_idx = compute_distance_to_clusters(x, centers, dist_type)

            if min_dist > max(radius[min_idx], dr):
                centers.append(x)
                NSamples.append(1.0)
                M.append(x*x)
                radius.append(0.0)
                idx_saver +=1
                idx.append(idx_saver)

            else:
                NSamples[min_idx] += 1
                centers[min_idx] = centers[min_idx]+(x-centers[min_idx])/(NSamples[min_idx])
                M[min_idx] = M[min_idx] + (x*x - M[min_idx])/(NSamples[min_idx])
                Sigma = M[min_idx] - centers[min_idx]*centers[min_idx]
                radius[min_idx] = np.sum(Sigma)
                idx.append(min_idx)


    ClusterStructureOut = {}
    ClusterStructureOut['centers'] = centers
    ClusterStructureOut['radius'] = radius
    ClusterStructureOut['NSamples'] = NSamples
    ClusterStructureOut['M'] = M
    ClusterStructureOut['idx'] = idx
    ClusterStructureOut['delta'] = delta



    if verbose:
        print("================= Results =====================")
        print("n_clusters={:d}".format(len(centers)))
        print("-- ClusterStructureOut['centers'] = ")
        print(centers)
        print("-- ClusterStructureOut['radius'] = ")
        print(radius)
        print("-- ClusterStructureOut['NSamples'] = ")
        print(NSamples)
        print("-- ClusterStructureOut['idx'] = ")
        print(idx)


    return ClusterStructureOut


