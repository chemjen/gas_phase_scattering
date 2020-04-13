import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from imageio import imread, mimsave
from tempfile import mkstemp, mkdtemp
tempdir = mkdtemp()
data = np.load('/reg/d/psdm/cxi/cxilo9215/results/jen/TMA_shorttime_3600batched2_data.npz')
x = data['q_bins']
y = data['phi_bins']
t = data['time_bins']*1000 + 1
Z = data['binned_pdiffs']
#x = np.delete(x, [0,28,48,49])
#Z = np.delete(Z, [0,28,48,49], axis=1)

X, Y = np.meshgrid(x, y)
cmap = cm.seismic
#cmap = cm.viridis
nbins = 31

#clim = max(abs(Z.min()),Z.max())
clim = 1
levels = MaxNLocator(nbins=nbins).tick_values(-1*clim, clim)
norm0 = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
pi=np.pi

frames = []

#print np.size(Z)
#print ((Z < 5)*(Z>-5)).sum()

for i in xrange(len(t)):
#	break
#	if i %5 !=0: continue
	z = Z[:,:,i]
	#figaspect = height/width
	fig = plt.figure(figsize = plt.figaspect(0.7))

	ax = fig.add_subplot(1,1,1)
	im = ax.pcolormesh(X, Y, z, cmap=cmap, norm=norm0)
	ax.set_xlabel('Time Delay (ps)', fontsize=20)
	ax.set_ylabel('$\phi$', fontsize=24)
	plt.yticks([-pi, -pi/2, 0, pi/2, pi],[r'-$\pi$', r'-$\frac{\pi}{2}$', '$0$', r'$\frac{\pi}{2}$', r'$\pi$'], size=20)
	ax.set_xlabel('q ($\AA^{-1}$)', fontsize=20)
	plt.xlim([x.min(), x.max()])
	ax.set_ylim([y.min(),y.max()])
	plt.yticks(size=18)
	plt.xticks(size=18)
	ax.set_title('%d fs delay' %t[i],fontsize=22)
	cb = fig.colorbar(im, ax=ax)
	cb.ax.tick_params(labelsize=14)
	cb.set_label('Percent Difference in Intensity', fontsize=18)
#	plt.show()
	_, filename = mkstemp(dir=tempdir)
	filename += '.png'
	plt.tight_layout()
	fig.savefig(filename)
	plt.close(fig)
	frames.append(filename)
#	del fig
#	if i > 10: break
	
images = []
for frame in frames:
    img = imread(frame)
    images.append(img)
mimsave('mygif.gif', images, duration = 0.25, loop=1)
#mimsave('temp.avi', images)
