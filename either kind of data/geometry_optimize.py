from pylab import *
from psana import *
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import minimize

#########################################################################################
# this script is to find the offsets in x0 and y0 (the center of the detector)
# and to find z0 (the detector distance)
# and to find phi0 (the offset in detector rotation)
# the I0 value is just a scaling factor

# Load premade .npy files
mask = np.load('/reg/d/psdm/cxi/cxilr8316/scratch/masks/run297_mask_xray.npy')
#x0 = (0 , 0, 8000, 0, 0.8)
x0 = np.array([  8.99930796e+02,  -1.79850420e+02,   8.46824843e+04, # shift 2
         7.36083925e-02,   1.0])
#bounds = None
bounds = ((-2000, 2000),(-2000, 2000),(60000, 100000), (0, np.pi), (0, 100))
method = 'SLSQP' #L-BFGS-B, SLSQP, TNC, BFGS (no bounds)
#method = 'TNC'
#method = 'L-BFGS-B'
#method = 'BFGS'

data = np.load('/reg/d/psdm/cxi/cxilr8316/scratch/sjen/CHD_longtime_3600batched_temp.npz')
image = data['offs']/data['offshots']
print data['offshots']
theory = np.loadtxt('/reg/d/psdm/cxi/cxilr8316/results/CHD_IAM.txt')

expimage = image*mask

ds = DataSource('exp=cxilr8316:run=120')
det = Detector('DscCsPad')
evt0 = ds.events().next()
xcoords = det.coords_x(evt0)  # coord in um
y = det.coords_y(evt0)
wavelength = 1.305 # xray wavelength, angstrom
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

def error(x):
	return (abs(experiment_image(x[0],x[1],x[2]) - theory_image(x[0],x[1],x[2],x[3],x[4]))**2).sum()
#calculate the absolute value of the error between the experimental and theoretical images

best_values = minimize(error, x0 = x0, method=method, bounds=bounds)
print method
print x0
print bounds
print best_values

#if best_values.success == True:
#    np.savetxt(‘best_values.txt’, best_values.x)

