#-------------------------------------------------------------------------------
from configuration import *
from scipy.integrate import odeint
import numpy as np
#-------------------------------------------------------------------------------
# see original http://stackoverflow.com/questions/19910189/scipy-odeint-with-complex-initial-values
def Odeintz(func, z0, t, **kwargs):
	"""An odeint-like function for complex valued differential equations."""
	# Disallow Jacobian-related arguments.
	_unsupported_odeint_args = ['Dfun', 'col_deriv', 'ml', 'mu']
	bad_args = [arg for arg in kwargs if arg in _unsupported_odeint_args]
	if len(bad_args) > 0:
		raise ValueError("The odeint argument %r is not supported by odeintz." % (bad_args[0],))

# Make sure z0 is a numpy array of type np.complex128.
	z0 = np.array(z0, dtype=np.complex128, ndmin=1)

	def realfunc(x, t, *args):
		z = x.view(np.complex128)
		dzdt = func(z, t, *args)
		# func might return a python list, so convert its return
		# value to an array with type np.complex128, and then return
		# a np.float64 view of that array.
		return np.asarray(dzdt, dtype=np.complex128).view(np.float64)

	result = odeint(realfunc, z0.view(np.float64), t, **kwargs)

	if kwargs.get('full_output', False):
		z = result[0].view(np.complex128)
		infodict = result[1]
		return z, infodict
	else:
		z = result.view(np.complex128)
		return z
#-------------------------------------------------------------------------------
