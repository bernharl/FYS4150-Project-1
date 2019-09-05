import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 3:
	raise RuntimeError("Only one data file can be plotted.\n Please enter only one data file in command line")


if len(sys.argv) == 3:
	filename_thomas = sys.argv[1]
	filename_special = sys.argv[2]

# if len(sys.argv) < 2:
# 	datafile = input("Please enter the name of the data file you wish to plot: ")



data_thomas = np.loadtxt(filename_thomas, skiprows = 1)
data_special = np.loadtxt(filename_special, skiprows = 1)
numerical_thomas = data_thomas[:, 0]
numerical_special = data_special[:, 0]
analytical = data_thomas[:, 1]
x = data_thomas[:, 2]


plt.plot(x, numerical_thomas, label = "Numerical, Thomas")
plt.plot(x, analytical, "r--", label = "Analytical")
plt.plot(x, numerical_special, "--", label = "Numerical, Special")
plt.legend()
plt.savefig("../doc/figures/graphs.eps")
