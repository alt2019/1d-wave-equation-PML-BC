# check python version
import sys
assert sys.version_info >= (3, 7)

# extra libs
import numpy as np

# needed project files 
from grid import Grid2D_TX
from wavequation import WaveEquation
from plotter import draw

L = 20

def phi(z, a = L / 5):
	return np.exp(-pow(z, 2) / (2 * pow(a, 2)))


if __name__ == '__main__':
	g = Grid2D_TX(
		x_start = -L, 
		x_end = L, 
		dx = 0.5/2, 
		t_start = 0.0, 
		t_end = 12.56*2*2, 
		dt = 0.005/2)

	we = WaveEquation(g, phi, c = 1)

	prec_sol = we.get_precise_solution().real
	draw(g, prec_sol[::50])	
	
	sol117 = we.solve_2nd_order_approximation_scheme()
	draw(g, sol117[::50])
	# abs_err = abs(sol117 - prec_sol)
	# plt.plot(g.X.real, np.transpose(abs_err[::100]))
	# plt.show()

	sol119 = we.solve_chess_scheme()
	draw(g, sol119[::50])
	# abs_err = abs(sol119 - prec_sol)
	# plt.plot(g.T[::100], abs_err[::100])
	# plt.plot(g.X.real, np.transpose(abs_err[::100]))
	# plt.show()

	sol121 = we.solve_PML(a = 1.0, delta = 2.10001)
	draw(g, sol121[::50])
	# abs_err = abs(sol121 - prec_sol)
	# plt.plot(g.X.real, np.transpose(abs_err[::100]))
	# plt.show()