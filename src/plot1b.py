import numpy as np
import matplotlib.pyplot as plt
import sys

""" Script loads and plots analytical and numerical solutions
	as function of x."""

if len(sys.argv) != 3:
	raise RuntimeError("Only one data file can be plotted.\n Please enter only one data file in command line")

if len(sys.argv) == 3:
	filename_thomas = sys.argv[1]
	filename_special = sys.argv[2]

# if len(sys.argv) < 2:
# 	datafile = input("Please enter the name of the data file you wish to plot: ")

"""Loading data"""
data_thomas_10 = np.loadtxt(filename_thomas % 10, skiprows = 1)
data_thomas_100 = np.loadtxt(filename_thomas % 100, skiprows = 1)
data_thomas_1000 = np.loadtxt(filename_thomas % 1000, skiprows = 1)
data_special_10 = np.loadtxt(filename_special % 10, skiprows = 1)
data_special_100 = np.loadtxt(filename_special % 100, skiprows = 1)
data_special_1000 = np.loadtxt(filename_special % 1000, skiprows = 1)

numerical_thomas_10 = data_thomas_10[:, 0]	# Numerical solution computed by Thomas algo. with n = 10
numerical_thomas_1000 = data_thomas_1000[:, 0]	# Numerical solution computed by Thomas algo. with n = 1000
numerical_special_10 = data_special_10[:, 0]	# Numerical solution computed by Specialised Thomas algo. with n = 10
numerical_special_1000 = data_special_1000[:, 0]	# Numerical solution computed by Specialised Thomas algo. with n = 1000
analytical = data_thomas_1000[:, 1]			# Analytical solution
x_10 = data_thomas_10[:, 2]					# x-values for n = 10
x_1000 = data_thomas_1000[:, 2]				# x-values for n = 1000

"""Ploting data"""
plt.plot(x_10, numerical_thomas_10, label = "Numerical, Thomas, n = 10")
plt.plot(x_1000, numerical_thomas_1000, label = "Numerical, Thomas, n = 1000")
plt.plot(x_1000, analytical, "--", label = "Analytical")
plt.plot(x_10, numerical_special_10, "--", label = "Numerical, Special, n = 10")
plt.plot(x_1000, numerical_special_1000, linestyle = "dashdot", label = "Numerical, Special, n = 1000")
plt.xlabel(r"$x$")
plt.ylabel(r"$u(x)$ $[v(x)]$")
plt.xticks(fontsize = 12)
plt.grid()
plt.legend(loc = 0)
plt.savefig("../doc/figures/graphs.eps")
