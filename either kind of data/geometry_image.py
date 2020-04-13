from pylab import *
from psana import *
from scipy.interpolate import InterpolatedUnivariateSpline
from radavg import RadialAverager 

#########################################################################################
# this script is to find the offsets in x0 and y0 (the center of the detector)
# and to find z0 (the detector distance)
# and to find phi0 (the offset in detector rotation)
# the I0 value is just a scaling factor
mask = np.load('/reg/d/psdm/cxi/cxilr8316/scratch/masks/run297_mask_xray.npy')
best_values = np.array([  8.93237364e+02,  -2.24028423e+02,   8.62315972e+04,
         9.06806988e-02,   8.35487644e-01])

data = np.load('/reg/d/psdm/cxi/cxilr8316/scratch/jen/CHD_longtime_3600batched_temp.npz')
image = data['offs']/data['offshots']
theory = np.loadtxt('/reg/d/psdm/cxi/cxilr8316/results/CHD_IAM.txt')
expimage = image*mask

ds = DataSource('exp=cxilr8316:run=120')
det = Detector('DscCsPad')
evt0 = ds.events().next()
xcoords = det.coords_x(evt0)  # coord in um
y = det.coords_y(evt0)
wavelength = 1.305 #wavelength, in Angstrom
 #Generate the theory image
ttheta = np.nan_to_num(np.arcsin((theory[:,0]*wavelength)/(4*np.pi))) #converting the q in the theory data to theta
tI = np.nan_to_num(theory[:,1]) #get rid of Not-a-Number values by converting them to zeros
tI = tI[np.nonzero(ttheta)] #get rid of the resulting zeros
ttheta = ttheta[np.nonzero(ttheta)] #make the x's and y's match
tspline = InterpolatedUnivariateSpline(ttheta,tI) #make a spline of the data
del tI, ttheta, theory #deletes variable that won't be used again, kind of unecessary but can be a good move when memory is an issue

# Define normalization
def norm(x): return np.nan_to_num(x/x.max()) #normalization/scaling scheme so that the experiment and theory can alig 

def theta(x0,y0,z0):
	return 0.5*np.arctan(sqrt((xcoords+x0)**2 + (y+y0)**2)/z0) # tan(2*theta) = r/z

def phi(x0, y0, phi0):
	return np.arctan2((y+y0),(xcoords+x0)) + phi0 	# phi = arctan(y/x)

def experiment_image(x0,y0,z0):
	return norm(expimage /np.cos(2*theta(x0,y0,z0))**3) 
#there's a cosine^2 factor for the dependence of scattering intensity on detector distance, 
# and another cosine factor for the effective area of the pixel

def theory_image(x0,y0,z0,phi0,I0):
	return norm(tspline(theta(x0,y0,z0)))*I0*(np.sin(phi(x0, y0, phi0))**2 + np.cos(phi(x0, y0, phi0))**2*np.cos(2*theta(x0,y0,z0))**2)*mask
# converting the unpolarized theoretical scattering pattern, to a scattering pattern for polarized x-rays

def the(x):
	return (abs(experiment_image(x[0],x[1],x[2]) - theory_image(x[0],x[1],x[2],x[3],x[4]))**2).sum()
#calculate the absolute value of the error between the experimental and theoretical images

exp = experiment_image(*best_values[:3])
theo = theory_image(*best_values)	

thetavals = theta(*best_values[:3])
ra = RadialAverager(theta(best_values[0], best_values[1], best_values[2]), mask, n_bins=37)
theta_bins = ra.bin_centers

subplot(1,2,1)
plot(theta_bins, ra(exp), label='experiment')
plot(theta_bins, ra(theo), label='theory')
legend()
xlabel('theta')
ylabel('intensity')

subplot(1,2,2)
plot(theta_bins, ra(exp-theo), label='residuals')
title('residuals: total = %.2f' %error(best_values))
xlabel('theta')
ylabel('intensity')
show()


