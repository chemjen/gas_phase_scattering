import matplotlib.pyplot as plt
import numpy as np
from psana import *
from IB import ImageBinner
from radavg import RadialAverager
from angavg import AngularAverager

## Before running this script, you need to do geometry optimization and error analysis
## Variables to change:
expdir = '/reg/d/psdm/cxi/cxilr8316/'
x0, y0, Z, phi0 = 8.93237364e+02,  -2.24028423e+02,   8.62315972e+04, 9.06806988e-02
#x0, y0, Z, phi0 =  8.97518584e+02,  -2.15728513e+02,   8.63635640e+04, 9.06932863e-02
n_phibins=48
n_qbins = 50
mask0 = np.load(expdir+'results/jen/run297_mask.npy')
wavelength = 1.305 #angstrom
time_bins = np.hstack((np.arange(-0.85, -0.05, 0.05), np.arange(-0.05, 0.92, 0.025)))
data = np.load('/reg/d/psdm/cxi/cxilr8316/results/jen/CHD_shorttime_smaller_tbins.npz')

## setting up variables
ds = DataSource('exp=cxilr8316:run=281')
det = Detector('DetInfo(CxiDs1.0:Cspad.0)', ds.env())
evt0 = ds.events().next()
X = det.coords_x(evt0) + x0	# coord in um
Y = det.coords_y(evt0) + y0
rr = np.sqrt((X)**2 + (Y)**2) 
theta = np.arctan2(rr, Z) / 2.0
q = 4. * np.pi * np.sin(theta) / wavelength 
phi_values = np.arctan2((Y), (X)) + phi0 
phi_values[phi_values < -1.*np.pi] += 2.*np.pi
ons = data['ons']
offs = data['offs']
onshots = data['onshots'] 
offshots = data['offshots']
## photon counting
Nons = np.copy(ons)
Noffs = np.copy(offs)

#averaging
for i in xrange(len(time_bins)):
	ons[:,:,:,i] *= mask0
	ons[:,:,:,i] /= onshots[i]
offs *= mask0
offs /= offshots
## error = sqrt(number of photons) = (summed image) / (ADU per photon)
## error_pdiff**2 = (dpdiff/dIon)**2 * (error_Ion)**2 + (dpdiff/dIoff)**2 * (error_Ioff)**2
## pdiff = (Ion - Ioff) * 100/ Ioff
## dpdiff/dIon = 100/Ioff
## dpdiff/dIoff = -100*Ion/Ioff**2
## --> error_pdiff = 100/|Ioff| * sqrt(error_Ion**2 + Ion**2*error_Ioff**2/Ioff**2)

## ib averages the images into bins, resulting in a (phi_bins x q_bins) array
ib = ImageBinner(q, phi_values, mask=mask0, q_bins=n_qbins, phi_bins=n_phibins)
pixel_counts = ib.pixel_counts
#np.save('pixel_counts.npy', pixel_counts)
ibpdiffs = np.zeros((n_phibins, n_qbins, len(time_bins)))
ibons = np.zeros_like(ibpdiffs)
iboffs = ib(offs)
print iboffs.sum()
for i in range(len(time_bins)):
    y1 = ib(ons[:,:,:,i])
    ibons[:,:,i] += y1
    y = np.nan_to_num((y1 - iboffs) * 100./iboffs)
    ibpdiffs[:,:,i] += y

error_offs = np.sqrt(ib(Noffs)/ib.pixel_counts)/offshots
error_ons = np.zeros((n_phibins, n_qbins, len(time_bins)))
for i in xrange(len(time_bins)):
    Nph = ib(Nons[:,:,:,i])/ib.pixel_counts
    error_ons[:,:,i] += np.sqrt(Nph)/onshots[i]

error_ibpdiffs = np.zeros_like(ibpdiffs)
A = 100./np.abs(iboffs)
B = error_offs**2/iboffs**2
for i in range(len(time_bins)):
    C = ibons[:,:,i]**2
    D = error_ons[:,:,i]**2
    error_ibpdiffs[:,:,i] += A * np.sqrt(B*C + D)

## Angular averaging, get a (phi_bins) array
aa = RadialAverager(phi_values, mask=mask0, n_bins=n_phibins)
aapdiffs = np.zeros((n_phibins, len(time_bins)))
aaons = np.zeros_like(aapdiffs)
aaoffs = aa(offs)
for i in range(len(time_bins)):
    y1 = aa(ons[:,:,:,i])
    aaons[:,i] += y1
    y = np.nan_to_num((y1 - aaoffs) * 100./aaoffs)
    aapdiffs[:,i] += y

error_offs = np.sqrt(aa(Noffs)/aa.pixel_counts)/offshots
error_ons = np.zeros((n_phibins, len(time_bins)))
for i in xrange(len(time_bins)):
    Nph = aa(Nons[:,:,:,i])/aa.pixel_counts
    error_ons[:,i] = np.sqrt(Nph)/onshots[i]

error_aapdiffs = np.zeros_like(aapdiffs)
A = 100./np.abs(aaoffs)
B = error_offs**2/aaoffs**2
for i in range(len(time_bins)):
    C = aaons[:,i]**2
    D = error_ons[:,i]**2
    error_aapdiffs[:,i] += A * np.sqrt(B*C + D)

## Radial Averaging, returns (q_bins) array
ra = RadialAverager(q, mask=mask0, n_bins=n_qbins)
rapdiffs = np.zeros((n_qbins, len(time_bins)))
raons = np.zeros_like(rapdiffs)
raoffs = ra(offs)
for i in range(len(time_bins)):
    y1 = ra(ons[:,:,:,i])
    raons[:,i] += y1
    y = np.nan_to_num((y1 - raoffs) * 100./raoffs)
    rapdiffs[:,i] += y

error_offs = np.sqrt(ra(Noffs)/ra.pixel_counts)/offshots #np.sqrt(ra(Noffs))/offshots
error_ons = np.zeros((n_qbins, len(time_bins)))
for i in xrange(len(time_bins)):
    Nph = ra(Nons[:,:,:,i])/ra.pixel_counts
    error_ons[:,i] += np.sqrt(Nph)/onshots[i]

error_rapdiffs = np.zeros_like(rapdiffs)
A = 100./np.abs(raoffs)
B = error_offs**2/raoffs**2
for i in range(len(time_bins)):
    C = raons[:,i]**2
    D = error_ons[:,i]**2
    error_rapdiffs[:,i] += A * np.sqrt(B*C + D)


new_timebins = np.zeros_like(time_bins)
bin_edge = time_bins[0] - (time_bins[1] - time_bins[0])
binedges = np.hstack((bin_edge, time_bins))

for i in xrange(len(new_timebins)):
	new_timebins[i] += ((binedges[i] + binedges[i+1])/2 )
print new_timebins
np.savez(expdir+'results/jen/CHD_shorttime_smaller_tbins_final.npz', \
	q_bins=ra.bin_centers, phi_bins=aa.bin_centers, time_bins=new_timebins, \
	ibpdiffs=ibpdiffs, errors=error_ibpdiffs, pixel_counts=pixel_counts, \
	rapdiffs=rapdiffs, error_rapdiffs=error_rapdiffs, offshots = offshots,\
	aapdiffs=aapdiffs, error_aapdiffs=error_aapdiffs, shots_per_bin=onshots)

import scipy.io as sio
sio.savemat('CHD_shorttime_smaller_tbins.mat', \
	{'q_bins': ra.bin_centers, 'phi_bins': aa.bin_centers, 'time_bins': new_timebins, \
	'percent_difference_intensities': ibpdiffs, 'errors': error_ibpdiffs, \
	'radial_pdiffs': rapdiffs, 'radial_pdiff_errors': error_rapdiffs, \
	'angular_pdiffs': aapdiffs, 'angular_pdiff_errors': error_aapdiffs, \
	'radial_uvons': raons, 'radial_uvoffs': raoffs, 'pixel_counts': pixel_counts, \
	'radial_uvon_errors': error_ons, 'radial_uvoff_errors': error_offs, \
	'uvoff_shots': offshots, 'shots_per_bin': onshots })


	

