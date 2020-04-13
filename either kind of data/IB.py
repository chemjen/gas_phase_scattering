
import numpy as np


"""
Example use:

>>> # assume image, radius_of_each_pixel, mask, are all arrays w/same shape
>>> ra = RadialAveragar(radius_of_each_pixel, mask, n_bins=101)
>>> radial_average = ra(image)
>>> ra.bin_centers # will give e.g. x-axis in a plot

"""

class ImageBinner(object):

	def __init__(self, q_values, phi_values, mask, q_bins=30, phi_bins=30):
		"""
		Parameters
		----------
		q_values : np.ndarray (float)
		    For each pixel, this is the momentum transfer value of that pixel
		phi_values : np.ndarray (float)
			For each pixel, the azimuthal angle w.r.t. xray polarization of that pixel
		k : np.ndarray (int)
    		A boolean (int) saying if each pixel is masked or not
		q_bins : int
    		The number of bins in q to employ. If `None` guesses a good value.
		phi_bins: int
			The number of bins in phi to employ.
		"""

		self.q_values = q_values
		self.phi_values = phi_values
		self.mask = mask
		self.q_bins = q_bins
		self.phi_bins = phi_bins

		self.q_range = self.q_values.max() - self.q_values.min()
		self.q_width = self.q_range / (float(q_bins))
		self.phi_range = self.phi_values.max() - self.phi_values.min()
		self.phi_width = self.phi_range / (float(phi_bins))
		self.n_bins = phi_bins * q_bins
			
		self._bin_assignments = np.zeros_like(q_values, dtype='int32')
		phi_assignments = np.floor( np.abs(self.phi_values - 1e-10 - self.phi_values.min())/self.phi_width).astype(np.int32)
		for i in xrange(phi_bins):
			phislice = self.q_values * (phi_assignments == i)
			temp_assignments = np.floor( np.abs(phislice - 1e-10 - q_values.min())/ self.q_width).astype(np.int32)
			bin_assignments_i = (temp_assignments + (i*q_bins)) * (phislice != 0)
			self._bin_assignments += bin_assignments_i
		self._normalization_array = (np.bincount( self._bin_assignments.flatten(), weights=self.mask.flatten() ) \
                                    + 1e-100).astype(np.float)
		assert self.n_bins >= self._bin_assignments.max() + 1, 'incorrect bin assignments'
		self._normalization_array = self._normalization_array[:self.n_bins]

		return
    
	def __call__(self, image):
		"""
		Bin pixel intensities by their momentum transfer.

		Parameters
		----------            
		image : np.ndarray
		    The intensity at each pixel, same shape as pixel_pos
		Returns
		-------
		bin_centers : ndarray, float
		    The q center of each bin.
		bin_values : ndarray, int
		    The average intensity in each bin. Each column represents a bin in q, each row a bin in phi.
		"""

		if not (image.shape == self.q_values.shape):
		    raise ValueError('`image` and `q_values` must have the same shape')
		if not (image.shape == self.mask.shape):
		    raise ValueError('`image` and `mask` must have the same shape')

		weights = image.flatten() * self.mask.flatten()
		bin_values = np.bincount(self._bin_assignments.flatten(), weights=weights)
		bin_values /= self._normalization_array
		
		#assert bin_values.shape[0] == self.n_bins
		if self.n_bins > len(bin_values):
			bin_values = np.hstack((bin_values, np.zeros(self.n_bins-len(bin_values))))
		
		return np.reshape(bin_values, (self.phi_bins, self.q_bins))
    
	@property
	def q_centers(self):
		return (np.arange(self.q_bins) + 0.5) * self.q_width + self.q_values.min()

	@property
	def phi_centers(self):
		return (np.arange(self.phi_bins) + 0.5) * self.phi_width + self.phi_values.min()

	@property
	def pixel_counts(self):
		arr = np.hstack((self._normalization_array, np.zeros(self.n_bins-len(self._normalization_array))))
		return np.reshape(arr, (self.phi_bins, self.q_bins)) 
