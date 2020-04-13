from psana import *
import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

### NEED TO UPDATE THE TIMETOOL PARAMETERS, MASK, KAPTON_AVG, and TIMEBINS
run_num = 265
lb = 16.3 # lower bound for adu per pixel
batch_size = 3600
diode_threshold = 22000.
diode_avg = 26000.
fee_min, fee_max = 9495.276, 9503.381
#pressure_min, pressure_max = 5.89, 6.2
#the real time bins are -5, 1, 3, 7, 10, 30, 70, 100, 300, 700, 1000, 3000, 7000]
time_bins = [-1, 2, 4, 8, 11, 31, 71, 101, 301, 701, 1001, 3001, 7001] #ps
datadir = '/reg/d/psdm/cxi/cxilr8316/scratch/jen/%dshots/' %batch_size

# Pull parameters from the data
ds = DataSource('exp=cxilr8316:run=%d'% run_num)
cspad = Detector('DscCsPad')
evr = Detector('evr1')
las_stg   = Detector('LAS:FS5:VIT:FS_TGT_TIME_DIAL') #ns
fee       = Detector('SIOC:SYS0:ML00:AO541')
#sample_pressure = Detector('CXI:MKS670:READINGGET') #pressure, in Torr
diode = Detector('Dg3Imp')
acqiris = Detector('Acqiris')
#sample_pressure = Detector('CXI:MKS670:READINGGET') #pressure, in Torr

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
		np.savez(datadir+'run%d_batch%d.npz' % (run_num, n/batch_size), 
		uv_offs=uv_offs, shots_per_bin=shots_per_bin, offshots=offshots, uv_ons=uv_ons)

for evt_num, evt in enumerate(ds.events()): #in python, a shot is called an "event" or "evt"
	if (evt_num%batch_size == 0) and (evt_num != 0):
		if on_shots.sum() > 0:
			savedata(evt_num)
		on_shots = np.zeros_like(on_shots)
		uvoffs = np.zeros_like(uvoffs)
		uvons = np.zeros_like(uvons)
		off_shots=0
	if evt_num%size != rank: continue

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
	time_delay = -1000.*las_stg(evt) 
	#evt_pressure = sample_pressure(evt)
	#if evt_pressure is None: continue
	#if (evt_pressure > pressure_max) or (evt_pressure < pressure_min): continue
	#print 'rawr'
	if uv_on :
		acqdata = np.asarray(acqiris(evt))
		uvmax = np.argmax(acqdata[0,3,:])
		uv_time = acqdata[1,3,uvmax]	
		if time_delay < 31: # time bins -5 to 30 ps
			if (uv_time < 4.35e-7) or (uv_time > 4.38e-7): continue
		elif time_delay < 101: # 70, 100 ps
			if (uv_time < 4.347e-7) or (uv_time > 4.38e-7): continue
		elif time_delay < 701: #300, 700 ps
			if (uv_time < 4.347e-7) or (uv_time > 4.376e-7): continue 
		elif time_delay < 1001: #1000 ps
			if (uv_time < 4.336e-7) or (uv_time > 4.375e-7): continue
		elif time_delay < 3001: #3000 ps
			if (uv_time < 4.325e-7) or (uv_time > 4.347e-7): continue
		else: #7000 ps
			if (uv_time > 4.31e-7): continue
		time_delay = np.array([-1000.*las_stg(evt)])	
		# Sort the shot into the correct time bin
		time_bin_index = int(np.digitize(time_delay,time_bins)) 
		#print 'rawr!'
		uvons[:,:,:,time_bin_index] += img 
		on_shots[time_bin_index] += 1 #uv-on counts in each time bin
	else: 
		uvoffs  += img
		off_shots += 1	

savedata((evt_num/batch_size +1)*batch_size)
MPI.Finalize()







