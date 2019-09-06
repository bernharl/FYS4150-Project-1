import numpy as np
import matplotlib.pyplot as plt
import sys

"""This script loads and plots the log10 of the maximum relative error
   between analytical and numerical solutions of ODE solved by 
   either Thomas or Specialised Thomas algorithm."""

if len(sys.argv) > 2:
	raise RuntimeError("Only one data file can be plotted.\n Please enter only one data file in command line")

if len(sys.argv) == 2:
	datafile = sys.argv[1]

if len(sys.argv) < 2:
	datafile = input("Please enter the name of the data file you wish to plot: ")

data = np.loadtxt(datafile, skiprows = 1)	# Loading data from file
n = data[:, 0]							  	# Grid sizes
error = data[:, 1]							# log10 of relative error

"""Plotting log10 of max relative error as function of grid size n""" 
plt.plot(np.log10(n), error)
plt.grid()
plt.xlabel(r"$\log_{10}{n}$", fontsize = 12)
plt.ylabel(r"$\log_{10}(\varepsilon_{max})$", fontsize = 12)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.savefig("../doc/figures/errorplot.eps")
