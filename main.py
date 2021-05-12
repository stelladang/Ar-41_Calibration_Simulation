import os
import numpy as np
import matplotlib.pyplot as plt
from Simulation import Simulation
if __name__ == "__main__":
	bins = np.arange(start=0, stop=6100, step=5)
	sim1 = Simulation()
	# print(sim1.start_cascade(verbose=True))
	E, fig, ax = sim1.plot_KE_graph(starting_energy = 6098.9, n_cascades = 1e4, title="Ar-41 Cascade KE Distribution", custom_bins = bins)
	ax.set_yscale("log")
	plt.show()

	### THERE'S AN OUTLIER AT 5582 keV with a density of 0.1168
	