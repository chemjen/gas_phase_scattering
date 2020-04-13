from psana import *
import numpy as np
import matplotlib.pyplot as plt
### NEED TO UPDATE THE TIMETOOL PARAMETERS, MASK, KAPTON_AVG, and TIMEBINS

run_num = 266
#265-280
#frontmask = np.load('/reg/d/psdm/cxi/cxilo9215/scratch/sjen/run84_mask_thr2.npy')
# Pull parameters from the data
ds = MPIDataSource('exp=cxilr8316:run=%d'% run_num)
front = Detector('DscCsPad') #front CSPAD
fee       = Detector('SIOC:SYS0:ML00:AO541') # FEE gas detector, give xray energy in eV
las_stg   = Detector('LAS:FS5:VIT:FS_TGT_TIME_DIAL') #ns
sample_pressure = Detector('CXI:MKS670:READINGGET') #pressure, in Torr
smldata = ds.small_data('/reg/d/psdm/cxi/cxilr8316/scratch/sjen/run%d_stats.h5' %run_num, gather_interval=1000)
evr = Detector('evr1')
diode = Detector('Dg3Imp')
acqiris = Detector('Acqiris')

# Reset variables for each batch of shots
for evt_num, evt in enumerate(ds.events()): #in python, a shot is called an "event" or "evt"
	#if evt_num > 10: break
#	print relative_time(evt)		
	# Sort shots according to EVR codes
	event_codes = evr(evt)
	if event_codes is None: continue
	uv_on = (183 in event_codes)
	xray_on = (162 not in event_codes)
	if not xray_on : continue
	img= np.copy(front.calib(evt)) #get the scattering pattern
	if img is None: continue
#	if img.sum() < detector_threshold: continue
	evt_diode = diode(evt)
	if evt_diode is None: continue
	diode1_intensity = evt_diode[1,:].max() - evt_diode[1,:].min()	
	diode2_intensity = evt_diode[3,:].max() - evt_diode[3,:].min()
	fee_data = fee(evt)
	if fee_data is None: continue
    # Get the time delay for the current shot
	time_delay = -1000.*las_stg(evt) 
	smldata.event(fee_data = fee_data, sample_pressure=sample_pressure(evt), diode1=diode1_intensity, diode2=diode2_intensity, img_intensity=img.sum())
	if uv_on :
		acqdata = np.asarray(acqiris(evt))
		uvmax = np.argmax(acqdata[0,3,:])
		uv_time = acqdata[1,3,uvmax]
	#	print evt_num, ':', uv_time
		# Sort the shot into the correct time bin
		smldata.event(uv_arrival=uv_time, relative_time=time_delay)
	else: continue 
		# Don't worry about time binning, just add the UV-off shot
smldata.save()







