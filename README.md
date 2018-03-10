# Temporal Human Action Segmentation via Dynamic Clustering

## Abstract
We present an effective and fast dynamic clustering algorithm for the task of temporal human action segmentation in videos, which has comprehensive applications such as robotics, motion analysis, patient monitoring and so forth. Our proposed algorithm is unsupervised, simple, and applicable in both the online and offline settings, where some standard post-processing method is applied to refine our result. To demonstrate the effectiveness, we perform extensive experiments of processing video streams from a generic multi-dimensional time series and different kinds of features. Our algorithm achieves the state-of-the-art results for both online and offline settings.

## Algorithm
Our dynamic clustering algorithm is a generic approach to processing sequential feature vectors (e.g. body pose from each frame) in an online manner. It consists of two steps: (1) initialization and (2) online updating.

__(1) Initialization__:
Initialization is completed via spectral clustering. In concrete, it captures the first batch of data and uses fully-connected graph spectra to determine the number of clusters in the current data batch. Then k-means is applied to determine initial cluster statistics.

__(2) Online update__:
Afterwards, our dynamic clustering algorithm processes each individual sample along the timestamps and update the cluster structure, which incorporates adding new clusters and updating cluster parameters. 

__Remarks__:
  * In contrast to many clustering algorithms, our method does not need to select k in advance. Instead, the cluster parameters and the number of clusters are derived jointly.

  * Other methods such as Dirichlet process mixture models (DPMMs) are also able to derive data-adaptive model complexity as well. However, DPMM allocates samples to clusters according to a Dirichlet process prior. Our method allocates samples to clusters depending on the data distance. Considering a video of 10 frames, the first 7 frames are
’walking’ and the following 2 frames are ’standing’. Due to the temporal coherence, the last frame tends to be
’standing’ as well. The DPMM model has larger probability to label this frame as ’walking’, since the ’walk-
ing’ cluster has more samples. But our algorithm will assign this frame to the ’standing’ cluster due to a lower
distance.

  * Our method shares several similarities with the adaptive resonance theory and hence works in a simplified manner of the brain.


## Experiments
Our dynamic clustering algorithm can play multiple roles to segment human actions temporally. 






Before running this matlab script, please ensure that:

  (1) Third-party libraries, such as TSC, KTC, ACA, matconvnet and vlfeat, have been installed.
  
  (2) The datasets have been downloaded and the interfaces have been setup. 
  
  (3) Video frame features from IDT+FV and VGG16 have been prepared.
  
  (4) The mex code is tested with Ubuntu 16, Matlab 2017a.

For any other inquiry, please contact
yan.zhang@uni-ulm.de
