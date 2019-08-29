import numpy as np
import matplotlib.pyplot as plt
import sys
	
if len(sys.argv) > 2:
	print("Only one data file can be plotted.\n Please enter only one data file in command line")
	sys.exit(1)
	
if len(sys.argv) == 2:
	datafile = sys.argv[1]
	
if len(sys.argv) < 2:
	datafile = input("Please enter the name of the data file you wish to plot: ")
	
	

data = np.loadtxt(datafile, skiprows = 1)
v = data[:, 0]
u = data[:, 1]
x = data[:, 2]
error = data[:, 3]

plt.plot(x, v, label = "Numerical")
plt.plot(x, u, "r--", label = "Analytical")
plt.legend()
plt.show()
