# check python version
import sys
assert sys.version_info >= (3, 7)

# extra libs
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib as mpl

# needed project files 
from grid import Grid2D_TX

def maximize_window(is_maximazed_window_needed: bool = False, 
					is_full_screen_needed: bool = False):
	"""
	This function maximizes window of plot 
	Parameters
	----------
	is_maximazed_window_needed : bool = False
		if True, the window of matplotlib figure will be maximized
	is_full_screen_needed: bool = False
		if True, the matplotlib figure will be full-screen
	"""
	manager = plt.get_current_fig_manager() 
	if is_maximazed_window_needed is True:
		manager.window.state('zoomed') # maximized window
	if is_full_screen_needed is True:
		manager.full_screen_toggle() # full-screen mode	


def draw(grid: Grid2D_TX, sol: np.ndarray, 
		 repeat_animation: bool = False,
		 is_maximazed_window_needed: bool = False, 
		 is_full_screen_needed: bool = False):
	"""
	This function implements time animation of solution 
	Parameters
	----------
	grid : Grid2D_TX
		grid with given function to animate
	sol : np.ndarray
		given function to animate
	repeat_animation: bool = False
		if True, animation will be done cyclically
	is_maximazed_window_needed : bool = False
		if True, the window of matplotlib figure will be maximized
	is_full_screen_needed: bool = False
		if True, the matplotlib figure will be full-screen
	"""
	X, dx, T, dt = grid.make_grid("real")

	plt.style.use('dark_background') 
	fig = plt.figure() 
	ax1 = fig.add_subplot(1,1,1, 
						  xlim = (X[0], X[-1]))

	def animate(t): 
		ax1.clear() 
		ax1.set_ylabel(r'$E(x)$') 
		ax1.set_xlabel(r'$X$') 
		ax1.set_ylim(0, 1.05)
		string =  't = {t}, T = {T:10f} '.format(
					 t = t,
					 T = (T[-1] - T[0]) / len(T) * t)
		text_pos_ID = int((len(X) - len(string)) / 2)
		ax1.text(X[text_pos_ID], 1.1, s = string)
		ax1.plot(X, abs(sol[t])) 
		ax1.grid(True) 

	anim = animation.FuncAnimation(
		fig = fig, 
		func = animate, 
		frames = len(sol), 
		interval = 1, 
		repeat = repeat_animation) 

	maximize_window(is_maximazed_window_needed, is_full_screen_needed)

	plt.show()


if __name__ == '__main__':
	... 