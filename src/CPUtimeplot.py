import numpy as np
import matplotlib.pyplot as plt
import sys

""" This script loads and plots the CPU time per algorithm as function of 
	grid size n."""

if len(sys.argv) > 3:
	print("Only two data files can be plotted.\n Please enter only two data files in command line")
	sys.exit(1)

if len(sys.argv) == 3:
	thomas_name = sys.argv[1]
	special_name = sys.argv[2] 

if len(sys.argv) < 2:
	datafile = input("Please enter the name of the data files you wish to plot: ")

"""Loading data"""
nLU = np.loadtxt("nLU.dat")
timeLU = np.loadtxt("timeLU.dat")
thomas_data = np.loadtxt(thomas_name, skiprows = 1)
special_data = np.loadtxt(special_name, skiprows = 1)
n = thomas_data[:, 0]					# Grid size
CPU_time_thomas = thomas_data[:, 2]		# CPU time of Thomas algo.
CPU_time_special = special_data[:, 2]	# CPU time of Specialized algo.

"""Ploting CPU times vs grid size n"""
plt.plot(np.log10(nLU), np.log10(timeLU), "ro" ,label = "LU decomposition")
plt.plot(np.log10(n), np.log10(CPU_time_thomas), label = "Thomas alg.")
plt.plot(np.log10(n), np.log10(CPU_time_special), label = "Specialized alg.")
plt.legend(loc = 0)
plt.xlabel(r"$\log_{10} n$")
plt.ylabel(r"$\log_{10} t_{CPU}$ [ms]")
plt.grid()
plt.savefig("../doc/figures/CPUtimeplot.eps")
