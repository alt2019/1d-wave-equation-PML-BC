# check python version
import sys
assert sys.version_info >= (3, 7)

# pre-installed libs
import math as m

# extra libs
import numpy as np

# needed project files 
from grid import Grid2D_TX

class WaveEquation:
	"""
	Class implements solution of 1d wave equation
	"""
	def __init__(self, 
				 grid: Grid2D_TX, 
				 init_func: np.ndarray,
				 c: float = 1.0):
		"""
		Parameters
		----------
		grid: Grid2D_TX
			Class grid to be used later
		init_func: np.ndarray
			Initial condition at zero time
		c: float
			Velocity of wave. Default is c = 1
		"""
		self.grid = grid
		self.init_func = init_func
		self.c = c

		# initialise grid
		self._set_grid()	
		# Courant number
		self.r = self.c * self.grid.dt / self.grid.dx
		if self.r >= 1: 
			print("Warning! Courant number c*dt/dx = {} > 1.\
				Scheme is unstable".format(self.r))
		else:
			print("Courant number c*dt/dx = {}".format(self.r))
		# initialize precise solution
		self._set_precise_solution()

	def _set_grid(self):
		self.X, self.dx, self.T, self.dt = self.grid.make_grid("real")

	def _set_precise_solution(self):
		"""
		This function creates 2d array of presice solution of given problem
		"""
		self.prec_sol = np.zeros((len(self.T), len(self.X)), dtype = complex)
		self.prec_sol[0] = self.init_func(self.X[:])
		for t in range(1, len(self.T)):
			self.prec_sol[t] = (self.init_func(self.X[:] - self.c * self.T[t-1])\
						    + self.init_func(self.X[:] + self.c * self.T[t-1])) / 2

	def get_precise_solution(self): return self.prec_sol

	def _sigma(self, z, a, delta):
		"""
		This function defines conductivity of free electric\
		 and magnetic charges (last mentioned are not observed\
		  in the real world) in PML scheme, where it is used\
		   as a non-reflecting boundary 
		Parameters
		----------
		a: float = 1.0
			attenuation parameter
		delta: float = 2.10001
			lenght of layer, where the solution of wave equation attenuate quickly
		"""		
		assert delta > 0
		if self.X[0] < z < self.X[0] + delta:
			return a * pow((z + abs(self.X[0]) - delta) / delta, 2)
		if self.X[0] + delta < z < self.X[-1] - delta:
			return 0
		if self.X[-1] - delta < z < self.X[-1]:
			return a * pow((z - self.X[-1] + delta) / delta, 2)

	def solve_2nd_order_approximation_scheme(self):
		"""
		This function implements solution of wave equation\
		 with 2nd order approximation scheme  
		"""	
		sol = np.zeros((len(self.T), len(self.X)))
		sol[0] = self.init_func(self.X[:])
		sol[1, 0] = 0; sol[1, -1] = 0
		sol[1, 1:-1:1] = sol[0, 1:-1:1]\
						 + pow(self.r, 2) / 2 *\
						  (sol[0, 2::1] - 2 * sol[0, 1:-1:1] + sol[0, 0:-2:1])
		for t in range(1, len(self.T)-1):
			sol[t+1, 0] = 0; sol[t+1, -1] = 0
			sol[t+1, 1:-1:1] = 2 * sol[t, 1:-1:1] - sol[t-1, 1:-1:1]\
				 		 + pow(self.r, 2) *\
				 		  (sol[t, 2::1] - 2 * sol[t, 1:-1:1] + sol[t, 0:-2:1])
		return sol

	def solve_chess_scheme(self):
		"""
		The main idea of th escheme is to replace wave equation\
		 by a system of Maxwell equations. The boundary conditions are\
		 such for ordinary Kauchy problem, so the solution must be the same as\
		 obtained in 2nd order approximation scheme   
		"""			
		E = np.zeros((len(self.T), len(self.X)))
		H = np.zeros((len(self.T), len(self.X)-1))
		E[0, :] = self.init_func(self.X[:])
		H[0, :] = 0.5 * self.r *\
			 (self.init_func(self.X[1::1]) - self.init_func(self.X[:-1:]))
		for t in range(0, len(self.T)-2):
			E[t+1, 0] = 0
			E[t+1, 1:-1:1] = E[t, 1:-1:1] + self.r * (H[t, 1::1] - H[t, 0:-1:1])
			E[t+1, -1] = 0
			H[t+1, 0::1] = H[t, 0::1] + self.r * (E[t+1, 1::1] - E[t+1, 0:-1:1])
		return E

	def solve_PML(self, a: float = 1.0, delta: float = 2.10001):
		"""
		This function implements algorythm for solving\
		 system of Maxwell equations with PML
		Parameters
		----------
		a: float = 1.0
		delta: float = 2.10001
			These parameters define stability\
			 of non-reflecting boundary conditions.\
			  They are transfered into self._sigma() function
		"""
		E = np.zeros((len(self.T), len(self.X)))
		H = np.zeros((len(self.T), len(self.X)-1))
		# creation of array sigma with double X-sturcture
		sigma = np.zeros(2 * len(self.X))
		g = self.grid
		G = Grid2D_TX(g.x_start, g.x_end, g.dx / 2)
		X2 = G.X.real
		for x in range(len(X2)):
			sigma[x] = self._sigma(X2[x], a = a, delta = delta)
		# solution
		E[0][:] = self.init_func(self.X[:])
		H[0, :] = 0.5 * self.r *\
			 (self.init_func(self.X[1:]) - self.init_func(self.X[:-1:]))
		for t in range(0, len(self.T)-2):
			E[t+1, 0] = 0
			E[t+1, 1:-2:1] = 1 / (1 + 2 * m.pi * g.dt * sigma[2:-4:2]) *\
								(E[t, 1:-2:1] * (1 - 2 * m.pi * g.dt * sigma[2:-4:2]) + 
									self.r * (H[t, 1:-1:1] - H[t, 0:-2:1]))
			E[t+1, -1] = 0
			H[t+1, 0:-2:1] = 1 / (1 + 2 * m.pi * g.dt * sigma[1:-5:2])  *\
						(H[t, 0:-2:1] * (1 - 2 * m.pi * g.dt * sigma[1:-5:2]) 
							+ self.r * (E[t+1, 1:-2:1] - E[t+1, 0:-3:1]))
		return E


if __name__ == '__main__':
	... 