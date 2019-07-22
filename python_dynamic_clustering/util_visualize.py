import numpy as np
import matplotlib.pyplot as plt

def show_seg(x, n_clusters):
	''' 
	This function is to show the segmentation, i.e. the cluster id sequence
	input:
		x 		    - the cluster id sequence, list
		n_clusters  - number of clusters, int
	'''
	x = list(x)
	n_frames = len(x)
	x_arr = np.array(x)
	x_color = x_arr/max(x)
	x_color = np.tile(x_color, [n_frames//10, 1])

	plt.figure()
	ax = plt.gca()

	ax.imshow(x_color, cmap = plt.get_cmap('viridis'))
	ax.set_xlabel('frames')
	ax.get_yaxis().set_visible(False)



