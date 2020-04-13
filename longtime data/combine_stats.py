import numpy as np
import h5py 
import matplotlib.pyplot as plt

runs = np.arange(265,281)

Fee = []
Diode1 = []
Diode2 = []
Pressure = []
Img_intensity = []
UV = []
Delay_time = []

for i in xrange(len(runs)):
	print runs[i]
	hf = h5py.File('/reg/d/psdm/cxi/cxilr8316/scratch/sjen/run%d_stats.h5' %runs[i], 'r')
	for item in hf:
		exec('%s = np.array(hf[item])' %str(item))
	Fee.extend(fee_data)
	Diode1.extend(diode1)
	Diode2.extend(diode2)
	Pressure.extend(sample_pressure)
	Img_intensity.extend(img_intensity)
	UV.extend(uv_arrival)
	Delay_time.extend(relative_time)

UV_arrival = np.nan_to_num(np.array(UV))
UV = UV_arrival[np.nonzero(UV_arrival)]
Delay_time = np.array(Delay_time)
Delay_time = Delay_time[np.nonzero(UV_arrival)]

def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
    """
    #arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.nanmedian(arr)
    return med - 3.*np.nanmedian(np.abs(arr - med)), med + 3.*np.nanmedian(np.abs(arr - med))

Fee = np.asarray(Fee)
print mad(Fee)
print np.median(Fee)
dd = Fee[Fee > mad(Fee)[0]]
ddd = dd[dd < mad(Fee)[1]]
print  np.mean(ddd), 'xray energy'
h =  4.135667516e-15 #Planck's constant, eV*s
c = 299792458. # m/s
print '=', h*c/np.mean(dd) * 1e10, 'angstrom'



plt.hist(Fee, 100)
plt.xlabel('xray energy (eV)')
plt.ylabel('counts')
plt.show()

plt.subplot(1,2,1)
plt.hist(Pressure)
plt.xlabel('sample pressure')
plt.ylabel('counts')
ax = plt.subplot(1,2,2)
plt.plot(Pressure)
plt.xlabel('shot')
ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
plt.ylabel('sample pressure')
plt.show()

ax = plt.subplot(2,2,1)
plt.hist(Diode1,100)
ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
plt.xlabel('diode 1 reading')
plt.ylabel('counts')
ax = plt.subplot(2,2,2)
plt.hist(Diode2,100)
ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
plt.xlabel('diode 2 reading')
plt.ylabel('counts')
ax = plt.subplot(2,2,3)
plt.scatter(Img_intensity, Diode1)
plt.ylabel('diode 1 reading')
ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
plt.xlabel('cspad intensity')
ax = plt.subplot(2,2,4)
plt.scatter(Img_intensity, Diode2)
plt.ylabel('diode 2 reading')
ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
plt.xlabel('cspad intensity')
plt.show()

delay_times = [-5, 1, 3, 7, 10, 30, 70, 100, 300, 700, 1000, 3000, 7000]
j = 1
for i in xrange(len(delay_times)):
	t = delay_times[i]
	indices = (Delay_time > (t - 1))*(Delay_time < (t + 1))
	uv_arrivals = UV[np.nonzero(indices)]
	plt.subplot(2,3,j)
	plt.hist(uv_arrivals, 100)
	plt.title('time = %d ps' %t)
	j += 1
	if j == 7:
		plt.tight_layout()	
		plt.show()
		j = 1
plt.show()
plt.hist(uv_arrivals, 100)
plt.title('time = %d ps' %t)
plt.show()
