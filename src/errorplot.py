import numpy as np
import matplotlib.pyplot as plt
import sys


if len(sys.argv) > 2:
	raise RuntimeError("Only one data file can be plotted.\n Please enter only one data file in command line")

if len(sys.argv) == 2:
	datafile = sys.argv[1]

if len(sys.argv) < 2:
	datafile = input("Please enter the name of the data file you wish to plot: ")

data = np.loadtxt(datafile, skiprows = 1)
n = data[:, 0]
error = data[:, 1]

plt.plot(np.log10(n), error)
plt.grid()
plt.xlabel(r"$\log_{10}{n}$")
plt.ylabel(r"$\varepsilon$")
plt.savefig("../doc/figures/errorplot.eps")
