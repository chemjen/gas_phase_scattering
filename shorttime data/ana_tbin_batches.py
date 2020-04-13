from psana import *
import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
### NEED TO UPDATE THE TIMETOOL PARAMETERS, MASK, KAPTON_AVG, and TIMEBINS

run_num = 288
if rank ==0:
	print run_num
lb = 16.3
batch_size = 3600
fee_min, fee_max = 9499.529, 9506.037
diode_threshold = 24000.
diode_avg = 27000.
TT_amp_min, TT_amp_max = 0.03, 0.14
TT_pos_min, TT_pos_max = 350., 650.
TT_fwhm_min, TT_fwhm_max = 40., 100.
datadir = '/reg/d/psdm/cxi/cxilr8316/scratch/jen/%dshots/' %batch_size
## TIMEBINNING
time_bins = np.hstack((np.arange(-1, -0.5, 0.05), np.arange(-0.5, 0.6, 0.025)))
timebinmin, timebinmax = time_bins.min(), time_bins.max()
timebinsize = time_bins[1] - time_bins[0]

# Pull parameters from the data
ds = DataSource('exp=cxilr8316:run=%d'% run_num)
cspad = Detector('DscCsPad')
evr = Detector('evr1')
las_stg = Detector('CXI:LAS:MMN:04.RBV') # laser stage in mm, time zero found at -20.6605 mm
tt_pos    = Detector('CXI:TTSPEC:FLTPOS')
tt_amp    = Detector('CXI:TTSPEC:AMPL')
tt_fwhm = Detector('CXI:TTSPEC:FLTPOSFWHM') # TT camera peak FWHM, in pixels 
fee       = Detector('SIOC:SYS0:ML00:AO541')
diode = Detector('Dg3Imp')
#sample_pressure = Detector('CXI:MKS670:READINGGET') #pressure, in Torr

# Time-tool calibration
def relative_time(evt): 
#    from docs >> fs_result = a + b*x + c*x^2, x is edge position
	a = -1.652557820602 
	b = 0.003833605540 #ps
	c = -0.000000740468
	x = tt_pos(evt)
	tt_correction = a + b*x + c*x**2
	time_delay = (las_stg(evt) - 72.695) * 2./0.3 # convert the mm stage position to ps (c =0.3 mm/ps)
	return -1.*(time_delay + tt_correction)

# Reset variables for each batch of shots
on_shots = np.zeros(len(time_bins))
off_shots = 0
shp = (32,185,388)
shpt = (32,185,388,len(time_bins))
uvons = np.zeros(shpt)
uvoffs = np.zeros(shp)

def savedata(n):
	shots_per_bin = np.empty_like(on_shots)
	uv_ons = np.empty_like(uvons)
	uv_offs = np.empty_like(uvoffs)
	comm.Reduce(on_shots, shots_per_bin)
	offshots = comm.reduce(off_shots)
	comm.Reduce(uvons, uv_ons)
	comm.Reduce(uvoffs, uv_offs)

	if rank==0:
		if shots_per_bin.sum() > 0.1:
			np.savez(datadir+'run%d_batch%d.npz' % (run_num, n/batch_size),\
				uv_offs=uv_offs, shots_per_bin=shots_per_bin, offshots=offshots, uv_ons=uv_ons)

for evt_num, evt in enumerate(ds.events()): #in python, a shot is called an "event" or "evt"
	#if evt_num/batch_size > 19.2: break
	#if evt_num > 8000: break
	if (evt_num%batch_size == 0) and (evt_num != 0):
		savedata(evt_num)
		on_shots = np.zeros_like(on_shots)
		uvoffs = np.zeros_like(uvoffs)
		uvons = np.zeros_like(uvons)
		off_shots=0

	if evt_num%size != rank: continue
	#if evt_num%200 == 0: print evt_num
	# Sort shots according to EVR codes
	event_codes = evr(evt)
	if event_codes is None: continue
	uv_on = (183 in event_codes)
	xray_on = (162 not in event_codes)
	if not xray_on : continue
	img= np.copy(cspad.calib(evt)) #get the scattering pattern
	if img is None: continue
#	if img.sum() < detector_threshold: continue
	fee_data = fee(evt)
	if (fee_data > fee_max) or (fee_data < fee_min): continue
	evt_diode = diode(evt)
	if evt_diode is None: continue
	diode2_intensity = evt_diode[3,:].max() - evt_diode[3,:].min()
	if diode2_intensity < diode_threshold: continue
#	evt_pressure = sample_pressure(evt)
#	if (evt_pressure < pressure_min) or (evt_pressure > pressure_max): continue
	try:
		img *= (img > lb)
	except TypeError:
		print('Error', evt_num, img, lb)
		continue
	img *= diode_avg/diode2_intensity
	time_delay = np.array([relative_time(evt)]) 
	if uv_on :
		if (tt_pos(evt) > TT_pos_max) or (tt_pos(evt) < TT_pos_min): continue 
		elif (tt_amp(evt) > TT_amp_max) or (tt_amp(evt) < TT_amp_min): continue
		elif (tt_fwhm(evt) > TT_fwhm_max) or (tt_fwhm(evt) < TT_fwhm_min): continue
		
		elif relative_time(evt) >= timebinmax: continue  #get rid of extras
		elif relative_time(evt) < (timebinmin - timebinsize): continue  
		# Sort the shot into the correct time bin
		time_bin_index = int(np.digitize(time_delay,time_bins)) 
		uvons[:,:,:,time_bin_index] += img 
		on_shots[time_bin_index] += 1 #uv-on counts in each time bin
#		smldata.event(tt_pos = tt_pos(evt), tt_amp=tt_amp(evt))
	else: 
		uvoffs  += img
		off_shots += 1	
		# Don't worry about time binning, just add the UV-off shot
#		continue
if evt_num%batch_size != 0:
	nbatch = (evt_num/batch_size +1)*batch_size
	savedata(nbatch)
MPI.Finalize()






