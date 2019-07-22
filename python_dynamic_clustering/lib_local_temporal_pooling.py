import numpy as np
import scipy
import scipy.ndimage.filters as filters
import scipy.signal as signal
from scipy.spatial.distance import cdist
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt

def soft_assignment(x, cluster_structure):

    '''
    feature encoding via soft assignment

    inputs:
        x - input sequence, list
        cluster_structure - output of dynamic clustering, dict

    return:
        output - encoded feature in array. Rows denote frames and columns denote feature dimensions
    '''

    centers = cluster_structure['centers']

    #compute pairwise distance between x & centers
    cluster_idx = cluster_structure['idx']
    centers_array = np.concatenate(centers,axis=0)
    dist_arr = cdist(x, centers_array, metric='sqeuclidean') 
    sim_arr = np.exp(-0.1*dist_arr)
    output = sim_arr/np.sum(sim_arr, axis=1)

    return output






def cal_motion_energy(x):
    '''
    motion energy = #transitions / len(x)

    input: 
        x - the local sequence, np array or list
    '''
    x = np.array(x)
    trans = np.sum(np.abs(np.diff(x))>1e-6)

    return float(trans)/x.shape[0]



def motion_energy_pooling(x, W, 
                          sigma = 1.0,
                          peak_distance=30,
                          plot_me=False):
    '''
    motion energy-based pooling: 
        (1) compute motion energy curve
        (2) detect peaks and valleys
        (3) find temporal boundaries, and assign 'moving actions' and 'still poses'

    input:
        x       - input sequence
        W       - time window to compute motion energy (an odd number)
        sigma   - for gaussian smoothing
        peak_distance - minimal distance between two peaks, we set to 1.5fps
    '''
    
    # compute the motion energy curve
    n_frames = x.shape[0]
    x_pad = np.pad(x, W//2, mode='edge')
    me_curve = [cal_motion_energy(x_pad[ii:ii+W]) for ii in range(n_frames)]
    me_curve = np.array(me_curve)

    # detect peaks and valleys
    ## first we remove noise using Gaussian smoothing
    ## note that this might remove small peaks
    me_curve = filters.gaussian_filter1d(me_curve, sigma=sigma)
    peak_idx, _ = signal.find_peaks(me_curve, distance=peak_distance)
    peak_width = signal.peak_widths(me_curve, 
                                    peak_idx, 
                                    rel_height=1.0)[0].astype(int)
    n_peaks = peak_idx.shape[0]

    
    # find temporal boundaries and assign labels
    
    action_bds = []
    for ii in range(n_peaks):
        action_bds.append([max(0,peak_idx[ii]-peak_width[ii]//2), 
                           min(peak_idx[ii]+peak_width[ii]//2, n_frames)
                           ] )

    action_labels = np.zeros(shape=x.shape)
    for peak_bd in action_bds:
        action_labels[peak_bd[0]:peak_bd[1]] = 1


    if plot_me:
        plt.figure()
        plt.plot(me_curve)

    return action_bds, action_labels





