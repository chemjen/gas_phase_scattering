from psana import *
import numpy as np
import matplotlib.pyplot as plt
### NEED TO UPDATE THE TIMETOOL PARAMETERS, MASK, and TIMEBINS

run_num = 281
#281-296
#frontmask = np.load('/reg/d/psdm/cxi/cxilo9215/scratch/sjen/run84_mask_thr2.npy')
#timebinmin, timebinmax = 5.5, 9.5
timebinmin, timebinmax = -1.2, 21
# Pull parameters from the data
ds = MPIDataSource('exp=cxilr8316:run=%d'% run_num)
front = Detector('DscCsPad') #front CSPAD
fee       = Detector('SIOC:SYS0:ML00:AO541') # FEE gas detector, give xray energy in eV
las_stg = Detector('CXI:LAS:MMN:04.RBV') # laser stage in mm, time zero found at -20.6605 mm
tt_pos    = Detector('CXI:TTSPEC:FLTPOS') # TT camera peak position, in pixels
tt_amp    = Detector('CXI:TTSPEC:AMPL') # TT camera peak amplitude, in pixels
tt_fwhm = Detector('CXI:TTSPEC:FLTPOSFWHM') # TT camera peak FWHM, in pixels 
sample_pressure = Detector('CXI:MKS670:READINGGET') #pressure, in Torr
smldata = ds.small_data('/reg/d/psdm/cxi/cxilr8316/scratch/sjen/run%d_stats.h5' %run_num, gather_interval=1000)
evr = Detector('evr1')
diode = Detector('Dg3Imp')

# Time-tool calibration
def relative_time(evt): 
#    from docs >> fs_result = a + b*x + c*x^2, x is edge position
# time_delay = (
	a = -1.652557820602 
	b = 0.003833605540 #ps
	c = -0.000000740468
	x = tt_pos(evt)
	tt_correction = a + b*x + c*x**2
	time_delay = (las_stg(evt) - 72.695) * 2./0.3 # convert the mm stage position to ps (c =0.3 mm/ps)
	return -1.*(time_delay + tt_correction)

# Reset variables for each batch of shots
for evt_num, evt in enumerate(ds.events()): #in python, a shot is called an "event" or "evt"
#	if evt_num > Nevents: break
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
	time_delay = np.array([relative_time(evt)]) 
	smldata.event(fee_data = fee_data, sample_pressure=sample_pressure(evt), diode1=diode1_intensity, diode2=diode2_intensity, img_intensity=img.sum())
	if uv_on :
		if relative_time(evt) > timebinmax: continue  #get rid of extras
		elif relative_time(evt) < timebinmin: continue  
		# Sort the shot into the correct time bin
		smldata.event(tt_pos = tt_pos(evt), tt_amp=tt_amp(evt), tt_fwhm = tt_fwhm(evt))
	else: continue 
		# Don't worry about time binning, just add the UV-off shot
smldata.save()







