import numpy as np
import scipy.io as sio
from lib_dynamic_clustering import *
from lib_local_temporal_pooling import *
import matplotlib.pyplot as plt
from util_visualize import *

'''
run hierarchical dynamic clustering for human motion segmentation:
Given an input feature sequence, we get segments boundaries of moving actions.
Dynamic clustering + motion energy-based pooling are employed.
'''


## read feature sequence
data_file_path = '/is/ps2/yzhang/Pictures/DynamicClustering/CMUMAD_features/CMUMAD_jointLocs_sub05_seq01.mat'
data = sio.loadmat(data_file_path)
X = data['X']

## run dynamic clustering
cluster_structure = dynamic_clustering(X, delta=1e-2, verbose=False)
cluster_idx = np.array(cluster_structure['idx'])


## run motion energy-based temporal pooling
action_bound, action_label = motion_energy_pooling(cluster_idx, 
                                                   W=30, 
                                                   sigma=2.5, 
                                                   peak_distance=30,
                                                   plot_me=True)

print(action_bound)

## visualize segmentation
# show_seg(cluster_structure['idx'], n_clusters=len(cluster_structure['centers']))
show_seg(action_label, n_clusters=2)
plt.show()


