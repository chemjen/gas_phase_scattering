from psana import *
import numpy as np
import matplotlib.pyplot as plt
### NEED TO UPDATE THE TIMETOOL PARAMETERS, MASK, KAPTON_AVG, and TIMEBINS

run_num = 280
fee_min, fee_max = 9495.276, 9503.381
diode_threshold = 22000
# Pull parameters from the data
ds = MPIDataSource('exp=cxilr8316:run=%d'% run_num)
fee       = Detector('SIOC:SYS0:ML00:AO541')
cspad = Detector('DscCsPad')
evr = Detector('evr1')
diode = Detector('Dg3Imp')
smldata = ds.small_data('/reg/d/psdm/cxi/cxilr8316/scratch/sjen/run%d_pixels.h5' %run_num, gather_interval=10000)

bins = np.linspace(-31,101,352)
counts = np.zeros(len(bins)-3)

# Reset variables for each batch of shots
for evt_num, evt in enumerate(ds.events()): #in python, a shot is called an "event" or "evt"
	# Sort shots according to EVR codes
	event_codes = evr(evt)
	if event_codes is None: continue
	uv_on = (183 in event_codes)
	xray_on = (162 not in event_codes)
	if not xray_on : continue
	img= np.copy(cspad.calib(evt)) #get the scattering pattern
	if img is None: continue
	evt_diode = diode(evt)
	if evt_diode is None: continue
	diode2_intensity = evt_diode[3,:].max() - evt_diode[3,:].min()
	if diode2_intensity < diode_threshold: continue
	fee_data = fee(evt)
	if fee_data is None: continue
	if (fee_data > fee_max) or (fee_data < fee_min): continue
	# Get the time delay for the current shot
	img[img >= 100] = 100.5
	img[img <= -30] = -30.5
	hist = np.histogram(img.flatten(),bins)
	counts += hist[0][1:350]

smldata.sum(counts)	
	# Don't worry about time binning, just add the UV-off shot
smldata.save(counts=counts)







