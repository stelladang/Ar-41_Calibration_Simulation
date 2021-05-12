import os
import numpy as np
import random
import operator
import matplotlib.pyplot as plt
from CurveFit import gaussian, fit_gauss

N_LINES = 4
SPEED_OF_LIGHT = 299792458.0 # in m/s

class EnergyLevel():

	def __init__ (self, file_name):

		start_level = 0.0 				# in keV
		half_life = 0.0					# in seconds
		transitions = np.array([])
		branching = np.array([])
		
		## opens a single data file and assigns values to instance fields
		with open(file_name) as file: 
			lines = [next(file).strip().replace(" ", "") for x in range(N_LINES)]

		for i in np.arange(N_LINES):
			parts = lines[i].split(':')
			if not (len(parts)==2):
				continue
			if(parts[0].lower()=="level"):
				self.start_level = float(parts[1])
			if(parts[0].lower()=="t12"):
				self.half_life = float(parts[1])
			if(parts[0].lower()=="transitions"):
				transitions_str = parts[1]
				self.transitions = [float(i) for i in transitions_str[1:len(transitions_str)-1].split(',')]
			if(parts[0].lower()=="branching"):
				branching_str = parts[1]
				self.branching = [float(i) for i in branching_str[1:len(branching_str)-1].split(',')]
	
	## for testing purposes - gives information about the energy level
	def to_string(self):
		print ("energy level: ", self.start_level, " keV")
		print ("half life: ", self.half_life, " seconds")
		print("transitions: ", self.transitions)
		print("branching ratios: ", self.branching)

class Simulation():

	def __init__(self, filepath=os.getcwd(), threshold=1e-10):
		
		self.hl_threshold = threshold
		self.energy_levels_dict = {}
		
		os.chdir(filepath) # not necessarily in order
		for file in os.listdir():
			filename, extension = os.path.splitext(file) 
			if extension == '.dat':
				curr_energy_level = EnergyLevel(file)
				self.energy_levels_dict[curr_energy_level.start_level] = curr_energy_level
	
	def start_cascade(self, cascade_start_energy=None, verbose = False): #-- pass in initial level and whether it is short lived

		if cascade_start_energy==None:
			cascade_start_level = self.energy_levels_dict[max(self.energy_levels_dict.keys())]
		else:
			cascade_start_level = self.energy_levels_dict[cascade_start_energy]

		curr_level = cascade_start_level

		tot_momentum = np.zeros(3)

		short_lived = True
		## uses a while loop to keep generating gammas while at a short-lived level
		while short_lived:
			## choose the next level using weighted probabilities
			next_level = self.energy_levels_dict[random.choices(population=curr_level.transitions, weights=curr_level.branching, k=1)[0]]
			gamma_energy = curr_level.start_level-next_level.start_level # in keV
			if verbose:
				print("gamma energy right now: ", gamma_energy, " keV")

			## calculate momentum using p = E/c
			p_mag = gamma_energy/SPEED_OF_LIGHT # in keV⋅s/m

			## make a random unit vector
			# -- pick two angles
			vec = np.array([random.gauss(0, 1) for i in range(3)])
			mag = sum(x**2 for x in vec) ** .5
			unit_vec = np.array([x/mag for x in vec])

			## multiply magnitude of momentum by unit vector
			p_vec = p_mag * unit_vec
			
			## add to total momentum vector
			tot_momentum += p_vec 
			tot_mag = (tot_momentum[0]**2 + tot_momentum[1]**2 + tot_momentum[2]**2) ** 0.5

			if verbose:
				print("total momentum right now: ", tot_momentum, "\nwith a magnitude of", tot_mag)
				print("therefore, the total KE right now should be ", tot_mag * SPEED_OF_LIGHT)
				print()

			if next_level.half_life >= self.hl_threshold:
				return tot_mag * SPEED_OF_LIGHT
				short_lived = False
			else:
				curr_level=next_level	

	def plot_KE_graph(self, starting_energy = None, n_cascades = 1e4, title="Total KE Distribution", custom_bins=None):
		KE_array = np.zeros(int(n_cascades))

		## run start_cascade() n_cascades times and collect result
		for i in np.arange(int(n_cascades)):
			KE_array[i]=self.start_cascade(cascade_start_energy=starting_energy)

		graph = plt.figure()
		freq, bins, patches = plt.hist(KE_array, density=True, bins=custom_bins)
		plt.title(title)
		plt.xlabel("Total KE [keV]")
		plt.ylabel("Density (normalized to unity)")
		# plt.yscale("log")
		bin_centers = custom_bins[1: len(custom_bins)]-(custom_bins[1]-custom_bins[0])/2
		# fit_gauss(bin_centers, freq, 50, 6100, fig=graph)
		plt.legend()
		ax0 = plt.gca()

		return KE_array, graph, ax0



