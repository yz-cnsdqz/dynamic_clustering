# Temporal Human Action Segmentation via Dynamic Clustering





## Software

Before running this matlab script, please ensure that:

  (1) Third-party libraries, such as TSC, KTC, ACA, matconvnet and vlfeat, have been installed.
  
  (2) The datasets have been downloaded and the interfaces have been setup. 
  
  (3) Video features __IDT+FV__ __VGG16__ __jointLocs__ __relativeAngle__ and __quaternions__ for all the used datasets (see below) have been prepared. One can download our prepared features via [this link](https://emotion.informatik.uni-ulm.de/public/DynamicClustering/). Perhaps you need to redefine the data path in the source files.
  
  (4) The mex code is tested with Ubuntu 16, Matlab 2017a. 
  
  ### NEW (5) In the python code, dynamic clustering and motion energy-based pooling are implemented.

## Paper and Citation
The manuscript can be found [here](https://arxiv.org/abs/1803.05790) and [here](http://bmvc2018.org/contents/papers/1000.pdf).

In case of using the code, please consider to cite following papers: 

       @article{1803.05790,
       Author = {Yan Zhang and He Sun and Siyu Tang and Heiko Neumann},
       Title = {Temporal Human Action Segmentation via Dynamic Clustering},
       Year = {2018},
       Journal = {arXiv preprint:1803.05790},
       }

      @conference{hdc:bmvc:2018,
        title = {Human Motion Parsing by Hierarchical Dynamic Clustering},
        author = {Zhang, Yan and Tang, Siyu and Sun, He and Neumann, Heiko},
        booktitle = {British Machine Vision Conference},
        month = sep,
        year = {2018},
        month_numeric = {9}
      }

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

  * Other methods such as Dirichlet process mixture models (DPMMs) are able to derive data-adaptive model complexity as well. However, DPMM allocates samples to clusters according to a Dirichlet process prior. Our method allocates samples to clusters depending on the data distance. Considering a video of 10 frames, the first 7 frames are
’walking’ and the following 2 frames are ’standing’. Due to the temporal coherence, the last frame tends to be
’standing’ as well. The DPMM model has larger probability to label this frame as ’walking’, since the ’walk-
ing’ cluster has more samples. But our algorithm will assign this frame to the ’standing’ cluster due to a lower
distance.

  * Our method shares several similarities with the adaptive resonance theory and hence works in a simplified manner of the brain.


## Experiments
Our dynamic clustering algorithm can play multiple roles to segment human actions temporally: (1) When the input is the frame-wise (or short-snippet-wise) feaure sequence, the output can be used as a codebook and then feature aggregation can be applied; (2)when the input is the action pattern sequence, its output is directly the set of action clusters. Therefore, we perform offline and online experiments separately, corresponding to the two aspects respectively.

### Datasets
(1) CMUMAD (http://www.humansensing.cs.cmu.edu/mad/)

(2) TUMKitchen (https://ias.in.tum.de/dokuwiki/software/kitchen-activity-data)

(3) MPII_HDM05_Sports (http://resources.mpi-inf.mpg.de/HDM05/)

### Features
(1) Improved dense trajectories + Fisher vectors

(2) the last layer of two-stream VGG16

(3) joint locations

(4) relative angles of joints

(4) quaternions of joints


### Comparison with state-of-the-art
* SC : spectral clustering
* TSC : temporal subspace clustering
* ACA : aligned cluster analysis
* EMS : efficient motion segmentation 
* DPMM : Dirichlet process mixture model
* DPMM-A : DPMM + temporal pooling + k-means
* ours : dynamic clustering + temporal pooling + k-means


(1) CMUMAD (precision/recall/runtime(in seconds))

| Algorithm | IDT+FV | VGG16 | JointLocation | RelativeAngle | Quaternion |
|-----------|-----------|-----------|-----------|-----------|-----------|
|SC| 0.57/0.85/203.3| 0.004/0.05/118.8 |0.02/0.13/113.0| 0.003/0.06/113.4 |0.01/0.11/125.3|
|TSC| 0.63/0.82/132.4| 0.01/0.2/38.2 |0.1/0.3/48.5| 0.05/0.29/41.6| 0.05/0.29/38.7|
|ACA|  __0.91__/0.83/547.7| 0.56/0.66/99.0| 0.55/0.68/221.5| 0.51/0.65/136.2| 0.55/__0.66__/168.8|
|EMS| 0.44/0.75/78.4| __0.67__/__0.73__/35.8| 0.34/0.78/33.0| 0.47/__0.89__/17.3| 0.6/0.51/9.2|
|DPMM| 0.4/0.73/507.8| 0.009/0.08/8.6| 0.02/0.12/17.8| 0.02/0.1/13.4| 0.02/0.11/11.6|
|DPMM-A| n/a |0.24/0.53/8.6| 0.37/0.54/17.8| 0.27/0.5/13.4| 0.39/0.58/11.6|
|ours| 0.56/__0.9__/__7.0__| 0.44/0.6/0.1| __0.82/0.86/0.1__| __0.63__/0.64/__0.1__| __0.63__/0.52/__0.1__|


![An example of segmentation is shown here.](https://github.com/IchBinYanZhang/dynamic_clustering/blob/master/CMUMAD_example_segmentation.png)

(2) TUMKitchen

details refer to manuscript

(3) MPII_HDM05_sports

details refer to manuscript



## Acknowledgement
We appreciate the kind technical supports from Viktor Kessler.
