import matplotlib.pyplot as plt
import numpy as np
import h5py
from psana import *
from radavg import RadialAverager

run = 316
shots = 30 #20000

ds = DataSource('exp=cxilr8316:run=28')
det = Detector('DetInfo(CxiDs1.0:Cspad.0)', ds.env())
evt0 = ds.events().next()
#mask0 = np.load('/reg/d/psdm/cxi/cxilr8316/scratch/masks/fbase_mask.npy')
mask0 = np.load('/reg/d/psdm/cxi/cxilr8316/scratch/masks/run157_mask_xray.npy')
nbins0 = 31
sf6 = False

filename = '/reg/d/psdm/cxi/cxilr8316/scratch/jen/run%d_%d_stats.h5' %(run, shots)
hf = h5py.File(filename, 'r')
for item in hf:
	print item
	exec('%s = np.array(hf[item])' %str(item))

print dark_shots, xray_shots

dark_front /= dark_shots
xray_front /= xray_shots

if sf6 is True:
	np.save('/reg/d/psdm/cxi/cxilr8316/scratch/run%d_sf6.npy' %run, xray_front)

X = det.coords_x(evt0)	# coord in um
Y = det.coords_y(evt0)
Xcenter = 0
Ycenter = 0
#Xcenter = (X.min() + X.max()) / 2.0
#Ycenter = (Y.min() + Y.max()) / 2.0

wavelength = 0.13051e-9	# X-ray wavelength, in meter
distance = 8e-2	# detector distance, in meter
rr = np.sqrt((X-Xcenter)**2 + (Y-Ycenter)**2) * 1e-6
theta = np.arctan2(rr, distance) / 2.0
q = 4. * np.pi * np.sin(theta) / wavelength * 1e-10

ra = RadialAverager(q, mask=mask0, n_bins=nbins0)
q_bins = ra.bin_centers[:-1]

### plot
def plot_radavg(data):
	plt.figure()
	plt.plot(q_bins, ra(data)[:-1])
	plt.show()

### 2D image
def plot_2dimage(data):
	img = det.image(evt0, data)
	plt.figure()
	plt.imshow(img, vmin=0, vmax=6)
	plt.colorbar()
	plt.show()

### update nbins
def update_nbins(nbins):
	global nbins0, ra, q_bins
	nbins0 = nbins
	ra = RadialAverager(q, mask=mask0, n_bins=nbins0)
	q_bins = ra.bin_centers[:-1]

### plot percent difference
def plot_pdiff(d1, d2):
	y1 = ra(d1)[:-1]
	y2 = ra(d2)[:-1]
	y = np.nan_to_num((y1-y2)*100/y2)
	plt.figure()
	plt.plot(q_bins, y)
	plt.show()

ax = plt.subplot(2,2,1)
plt.hist(diode1_intensity, 100)
plt.xlabel('diode intensity - Dg3Imp[1]')
plt.ylabel('counts')
plt.title('run %d' %run)
ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
ax = plt.subplot(2,2,2)
plt.hist(diode2_intensity, 100)
plt.xlabel('diode intensity - Dg3Imp[3]')
ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
plt.subplot(2,2,3)
plt.hist(sample_pressure, 20)
plt.xlabel('sample pressures')
ax = plt.subplot(2,2,4)
plt.hist(front_intensity,100)
plt.xlabel('CSPAD intensity')
ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
plt.show()

plt.hist(xray_energy,100)
plt.xlabel('xray energy (fee)')
plt.show()


plt.subplot(2,1,1)
plt.scatter(diode2_intensity, diode1_intensity)
plt.xlabel('diode - Dg3Imp[3]')
plt.ylabel('diode - Dg3Imp[1]')
plt.subplot(2,1,2)
plt.scatter(diode2_intensity, front_intensity)
plt.xlabel('diode - Dg3Imp[3]')
plt.ylabel('cspad intensity')
plt.show()





