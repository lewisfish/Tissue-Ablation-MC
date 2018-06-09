import matplotlib.pylab as plt
import numpy as np


x, y, z = np.loadtxt("depth.dat", unpack=True)
y *= 2.*.55/104.
plt.plot(x, y)
plt.xlabel("Power/W")
plt.ylabel("Depth/cm")

plt.show()
