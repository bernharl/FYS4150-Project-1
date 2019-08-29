import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data.dat", skiprows = 1)
v = data[:, 0]
u = data[:, 1]
x = data[:, 2]
error = data[:, 3]

plt.plot(x, v, label = "Numerical")
plt.plot(x, u, "r--", label = "Analytical")
plt.legend()
plt.show()
