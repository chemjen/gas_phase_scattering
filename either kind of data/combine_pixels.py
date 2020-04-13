import numpy as np
import h5py 
import matplotlib.pyplot as plt

runs = np.arange(265,281)
bins = np.linspace(-31,101,352)
Counts = np.zeros(349)
bin_centers = np.linspace((bins[0] + bins[1])/2, (bins[-2] + bins[-3])/2, 349)

for i in xrange(len(runs)):
	print runs[i]
	hf = h5py.File('/reg/d/psdm/cxi/cxilr8316/scratch/sjen/run%d_pixels.h5' %runs[i], 'r')
	for item in hf:
		exec('%s = np.array(hf[item])' %str(item))
	if i == 0:	
		Counts += counts

plt.plot(bin_centers, Counts)
plt.xlabel('ADU per pixel')
plt.ylabel('counts')
plt.show()
