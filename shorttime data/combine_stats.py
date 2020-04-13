import numpy as np
import h5py 
import matplotlib.pyplot as plt

#runs = np.arange(281,297)
runs = np.arange(282,297)

Fee = []
Diode1, Diode2, Pressure, Img_intensity = [], [], [], []
TT_pos, TT_amp, TT_fwhm = [], [], []
		
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
	TT_pos.extend(tt_pos)
	TT_amp.extend(tt_amp)
	TT_fwhm.extend(tt_fwhm)


TT_pos = np.array(TT_pos)
TT_amp = np.array(TT_amp)
TT_fwhm = np.array(TT_fwhm)
TT_pos = TT_pos[np.nonzero(np.nan_to_num(TT_pos))]
TT_amp = TT_amp[np.nonzero(np.nan_to_num(TT_amp))]
TT_fwhm = TT_fwhm[np.nonzero(np.nan_to_num(TT_fwhm))]

def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
    """
    #arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.nanmedian(arr)
    return med - 3.*np.nanmedian(np.abs(arr - med)), med + 3.*np.nanmedian(np.abs(arr - med))

print mad(Fee)

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
plt.tight_layout()
plt.show()

ax = plt.subplot(1,3,1)
plt.hist(TT_pos, 100)
plt.ylabel('counts')
ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
plt.xlabel('TT pixel position')
ax = plt.subplot(1,3,2)
plt.hist(TT_amp, 100)
plt.xlabel('TT peak amplitude')
ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
plt.subplot(1,3,3)
plt.hist(TT_fwhm, 100, range=[0, 150])
plt.xlabel('TT gaussian fwhm')
plt.tight_layout()
plt.show()
