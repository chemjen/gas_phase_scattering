from psana import *
import numpy as np

## Script for doing analysis of front and back detectors, and seeing the kapton ring, for building masks.

run = 279
Nevents =30000 #the maximimum number of shots you want to analyze
#ds = MPIDataSource('exp=cxilr8316:run=%d:dir=/reg/d/ffb/cxi/cxilr8316/xtc'% run) 
ds = MPIDataSource('exp=cxilr8316:run=%d' % run) #parallelized "small data" object, best to read the psana docs
front = Detector('DscCsPad', ds.env()) #the front detector
frontmask = np.load('/reg/d/psdm/cxi/cxilr8316/scratch/masks/fbase_mask.npy') # mask for front detector, can use a real mask once we make one
#lower_bound = -30. # lower bound for ADU per pixel
evr = Detector('evr1')
diode = Detector('Dg3Imp')
fee       = Detector('SIOC:SYS0:ML00:AO541')
sample_pressure = Detector('CXI:MKS670:READINGGET')
smldata = ds.small_data('/reg/d/psdm/cxi/cxilr8316/scratch/run%d_%d_stats.h5' %(run, Nevents/1000.), gather_interval=800) #prepare file to save

# EVR codes
LASERON		= 183
LASEROFF	= 184
XRAYOFF		= 162
XRAYOFF1	= 163

#Prepare variables, counting shots and saving images
xray_shots = 0
dark_shots = 0
xray_front = np.zeros_like(frontmask)
dark_front = np.zeros_like(frontmask)

#looping over the shots
for n, evt in enumerate(ds.events()):
#	if n > Nevents: break
	if evt is None: 
		print 'evt'
		continue
	img = front.calib_data(evt) # #grabbing detector images
	evt_diode = diode(evt)
	if evt_diode is None:
		print 'diode'
		continue
	diode1_intensity = evt_diode[1,:].max() - evt_diode[1,:].min()	
	diode2_intensity = evt_diode[3,:].max() - evt_diode[3,:].min()
	fee_data = fee(evt)
	evt_pressure = sample_pressure(evt)	
#ipm_data = IPM(evt)
	#sometimes a bad shot will be saved, without any images
	if img is None: 
		print 'img'
		continue
	img *= frontmask	
	evrcodes = evr(evt)
	if evrcodes is None:
		print 'evrcodes'
		continue
	dark = XRAYOFF in evrcodes
	if not dark:
		xray_front += img # * (img > lower_bound)
		xray_shots += 1
		smldata.event(front_intensity = img.sum(), diode1_intensity=diode1_intensity, diode2_intensity=diode2_intensity, sample_pressure = evt_pressure, xray_energy = fee_data)
#		smldata.event(fee_data = fee_data)
#		smldata.event(ipm_data = ipm_data)
	elif dark:
		dark_front += img  # * float(img > lower_bound)
		dark_shots += 1
# summing up the images and shots across ranks
xrayshots = smldata.sum(xray_shots)
darkshots = smldata.sum(dark_shots)
front_xrays = smldata.sum(xray_front)
front_dark = smldata.sum(dark_front)
# saving the summed values
smldata.save(xray_front=front_xrays, dark_front=front_dark, xray_shots=xrayshots, dark_shots=darkshots)

