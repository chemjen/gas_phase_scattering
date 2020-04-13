from masking import maskmaker
import numpy as np
import matplotlib.pyplot as plt

run = 262
Nevents = 20000
datadir = '/reg/d/psdm/cxi/cxilr8316/scratch/'
detector = 'front' #For LR83 we are only looking at the front detector
filename = datadir+'path/%dfile_%d.h5' %(run, Nevents/1000)
#whichmask = 'dark' #'dark' or 'xray'
whichmask = 'xray'

if whichmask == 'xray':
	oldmask = np.load(datadir+'masks/run%d_mask_dark.npy' %run)	
if whichmask == 'dark':
	oldmask = np.load(datadir+'masks/fbase_mask.npy') #mask of just edge & unbonded pixels

#maskmaker makes the mask, starting from the 'oldmask' input
mm = maskmaker('cxilr8316', run, detector, filename, oldmask)

xray_image = np.copy(mm.xray)
dark_image = np.copy(mm.dark)

#plt.hist(mm.fee)
#plt.title('histogram of xray energy')
#plt.show()
#plt.hist(mm.ipm)
#plt.show()

print 'average detector signal', np.average(mm.intensity)
plt.hist(mm.intensity)
plt.title('cspad intensity')
plt.ylabel('shot counts')
plt.xlabel('intensity per shot')
plt.show()

print 'use the histogram to find lower and upper bounds'
if whichmask == 'xray':
	plt.hist(xray_image.flatten(), bins=200, range = (-50, 150))
	plt.title('xray-on image')
if whichmask == 'dark':
	plt.hist(dark_image.flatten(), bins=100)
	plt.title('xray-off image')
plt.ylabel('pixel counts')
plt.xlabel('ADU per pixel')
plt.show()

print 'input the lower bound:'
lb = float(raw_input())
print 'input the upper bound:'
ub = float(raw_input())

assert ub > lb

name = datadir+'/masks/run%d_mask_%s.npy' %(run, whichmask)
newmask = mm(whichmask, lowerbound=lb, upperbound=ub, tolerance=2.5, name=name)
print newmask.sum()

img = mm.det.image(mm.evt0, xray_image)
plt.imshow(img, vmin=lb, vmax=ub)
plt.colorbar()
plt.title('xray-on cspad image')
plt.show()

img = mm.det.image(mm.evt0, dark_image)
plt.imshow(img, vmin=-2, vmax=2)
plt.title('dark cspad image')
plt.colorbar()
plt.show()

img = mm.det.image(mm.evt0, newmask* xray_image)
plt.imshow(img, vmin=lb, vmax=ub)
plt.title('dark cspad image')
plt.colorbar()
plt.title('masked cspad')
plt.show()
