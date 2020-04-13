import numpy as np

datadir = '/reg/d/psdm/cxi/cxilr8316/'
mask0 = np.load(datadir+'scratch/masks/run297_mask_xray.npy')

batch_size = 3600
time_bins = [-5, 1, 3, 7, 10, 30, 70 ,100, 300, 700, 1000, 3000, 7000]

offs = np.zeros((32,185,388))
offshots = 0
runs = np.arange(265, 281)
num_batches = 23517/batch_size + 1
print num_batches

for i in xrange(len(runs)):
	for batch in xrange(num_batches):
		try: 
			batch_data = np.load(datadir+'scratch/sjen/%dshots/run%d_batch%d.npz' \
				%(batch_size, runs[i], batch+1))
		except: 
			print runs[i], batch
			continue

		offs += batch_data['uv_offs']
		offshots += batch_data['offshots']

offs *= mask0/offshots
print offs.sum()


