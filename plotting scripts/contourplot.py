import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

data = np.load('/reg/d/psdm/cxi/cxilo9215/results/jen/TMA_shorttime_1800batched2_data.npz')
y = data['phi_bins']
x = data['time_bins']
z = data['angular_pdiffs']
#y = data['q_bins']
#x = data['radial_pdiffs']
#y = np.delete(y, [0,28,48,49])
#z = np.delete(z, [0,28,48,49], axis=0)

X, Y = np.meshgrid(x, y)

## This section sets up the colormap
## I have the color map set so it goes from -clim to +clim,
## diverging from 0 = white to blue = negative and red = positive
cmap = cm.seismic 
nbins = 41 #number of color bins
#clim = max(abs(z[:,:].min()),z[:,:].max())
clim = 0.5
levels = MaxNLocator(nbins=nbins).tick_values(-1.*clim, clim)
norm0 = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

pi=np.pi

fig = plt.figure(figsize = plt.figaspect(0.7))
ax = fig.add_subplot(1,1,1)
im = ax.pcolormesh(X, Y, z, cmap=cmap, norm=norm0)
ax.set_xlabel('Time Delay (ps)', fontsize=20)
ax.set_ylabel('$\phi$', fontsize=24)
plt.yticks([-pi, -pi/2, 0, pi/2, pi],[r'-$\pi$', r'-$\frac{\pi}{2}$', '$0$', r'$\frac{\pi}{2}$', r'$\pi$'], size=20)
#ax.set_ylabel('q ($\AA^{-1}$)', fontsize=20)
plt.xlim([-1, 5.3])
ax.set_ylim([y.min(),y.max()])
plt.yticks(size=18)
plt.xticks(size=18)
#ax.set_title('Radial Scattering',fontsize=22)
ax.set_title('Angular Scattering', fontsize=22)
cb = fig.colorbar(im, ax=ax)
cb.ax.tick_params(labelsize=14)
cb.set_label('Percent Difference in Intensity', fontsize=18)
plt.show()










