import matplotlib.pyplot as plt
import numpy as np
from psana import *
from radavg import RadialAverager
from angavg import AngularAverager
from IB import ImageBinner

ds = DataSource('exp=cxilr8316:run=12')
datadir = '/reg/d/psdm/cxi/cxilr8316/'
det = Detector('DetInfo(CxiDs1.0:Cspad.0)', ds.env())
evt0 = ds.events().next()
nbins0 = 50
mask0 = np.load(datadir+'scratch/masks/run297_mask_xray.npy')

batch_size = 3600
time_bins = [-5, 1, 3, 7, 10, 30, 70 ,100, 300, 700, 1000, 3000, 7000]

ons = np.zeros((32,185,388,len(time_bins)))
offs = np.zeros((32,185,388))
onshots = np.zeros(len(time_bins))
offshots = 0
runs = np.arange(265, 281)
uvoff_average = 1625167.46957
num_batches = 23517/batch_size + 1

for i in xrange(len(runs)):
	if runs[i] == 133: continue
	for batch in xrange(num_batches):
		try: 
			batch_data = np.load(datadir+'scratch/sjen/%dshots/run%d_batch%d.npz' \
				%(batch_size, runs[i], batch+1))
		except: 
			print runs[i], batch
			continue
			
		scaling_factor = uvoff_average * batch_data['offshots'] / (mask0 * batch_data['uv_offs']).sum()
		offs += batch_data['uv_offs']*scaling_factor
		ons += batch_data['uv_ons']*scaling_factor
		onshots += batch_data['shots_per_bin']
		offshots += batch_data['offshots']

np.savez(datadir+'scratch/sjen/CHD_longtime_%dbatched_temp.npz' %batch_size, \
		offs=offs, ons=ons, onshots=onshots, offshots=offshots)

