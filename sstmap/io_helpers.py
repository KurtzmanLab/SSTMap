# encoding: utf-8
import numpy as np

def are_you_numpy(a):

	"""
	Returns True if a is an instance of numpy.
	False otherwise.
	"""

	return type(a).__module__ == np.__name__


def make_grid(arrays, out=None):
	"""
	!!! Adapted from:
	!!! http://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations-of-two-arrays

	Generate a cartesian product of input arrays.

	Parameters
	----------
	arrays : list of array-like
		1-D arrays to form the cartesian product of.
	out : ndarray
		Array to place the cartesian product in.

	Returns
	-------
	out : ndarray
		2-D array of shape (M, len(arrays)) containing cartesian products
		formed of input arrays.

	Examples
	--------
	>>> make_grid(([1, 2, 3], [4, 5], [6, 7]))
	array([[1, 4, 6],
			[1, 4, 7],
			[1, 5, 6],
			[1, 5, 7],
			[2, 4, 6],
			[2, 4, 7],
			[2, 5, 6],
			[2, 5, 7],
			[3, 4, 6],
			[3, 4, 7],
			[3, 5, 6],
			[3, 5, 7]])

	"""

	arrays = [np.asarray(x) for x in arrays]

	dtype  = arrays[0].dtype

	n = np.prod([x.size for x in arrays])

	if out is None:

		out = np.zeros([n, len(arrays)], dtype=dtype)

	m = n / arrays[0].size

	out[:,0] = np.repeat(arrays[0], m)

	if arrays[1:]:

		make_grid(arrays[1:], out=out[0:m,1:])

		for j in xrange(1, arrays[0].size):

			out[j*m:(j+1)*m,1:] = out[0:m,1:]

	return out


def bounding_box_frac(frac_structure, delta=np.ones(3), _buffer=0., verbose=False):

	"""
	Input is structure in cart. or frac. coordinates as
	nx3 array (n= number of coordinates).
	Output is coordinate meshgrid array with coordinates of
	bounding box lattice as integers.
	"""

	bounding_min = np.array( [ np.min(frac_structure[:,0]),
							   np.min(frac_structure[:,1]),
							   np.min(frac_structure[:,2]) ], dtype=int )

	bounding_max = np.array( [ np.max(frac_structure[:,0]),
							   np.max(frac_structure[:,1]),
							   np.max(frac_structure[:,2]) ], dtype=int )

	bounding_min -= int(np.round(_buffer))
	bounding_max += int(np.round(_buffer))

	if verbose:
		print "Bounding min. ", bounding_min
		print "Bounding max. ", bounding_max
		print np.arange(bounding_min[2], bounding_max[2]+1, delta[2], dtype=int )

	return  make_grid ( [ np.arange(bounding_min[0], bounding_max[0]+1, delta[0], dtype=int ),
			  			  np.arange(bounding_min[1], bounding_max[1]+1, delta[1], dtype=int ),
			  			  np.arange(bounding_min[2], bounding_max[2]+1, delta[2], dtype=int ) ] )

