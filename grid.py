# check python version
import sys
assert sys.version_info >= (3, 7)

# extra libs
import numpy as np

class Grid2D_TX:
	"""
	Class implements grid over space and time 
	"""
	def __init__(self,
				 x_start: float = 0.0, x_end: float = 1.0, dx: float = 0.001,
				 t_start: float = 0.0, t_end: float = 1.0, dt: float = 0.001):
		"""
		Parameters
		----------
		x_start, x_end : float
			Left and right coordinates of cut over X-space
			default: x_start = 0.0, x_end = 1.0
		dx : float
			Step of grid over X-space cut
			default: dx = 0.001
		t_start, t_end : float
			Left and right coordinates of cut over T-space
			default: t_start = 0.0, t_end = 1.0
		dt : float
			Step of grid over T-space cut
			default: dt = 0.001
		"""
		self.x_start = x_start
		self.x_end = x_end
		self.dx = dx
		self.t_start = t_start
		self.t_end = t_end
		self.dt = dt

		# iniializing grid
		self._set_grid()

	def _set_grid(self):
		"""
		This function creates arrays over X, T in the form [start, stop, dt].
		X-array is complex, T-array is real
		"""
		self.X = np.arange(self.x_start, self.x_end + self.dx, self.dx, dtype = complex)
		self.T = np.arange(self.t_start, self.t_end + self.dt, self.dt)

	def make_grid(self, dtype = "complex"):
		"""
		This function returns arrays over space X - where X\
		 can be returned as complex or real, time T and its steps
		Parameters
		----------
		dtype : string
			Must be "complex"/"real"
		"""
		if dtype == "real":
			X =  self.X.real
		elif dtype == "complex":
			X = self.X
		elif dtype is not ("real" or "complex"):
			raise "dtype of returned arrays is not understood.\
				   Please specify 'comlex/real'"
		return X, self.dx, self.T, self.dt


if __name__ == '__main__':
	... 