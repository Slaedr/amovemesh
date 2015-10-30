#! /usr/bin/env python3
import sys
import numpy as np
from matplotlib import pyplot as plt

#fname = "points.dat"
fname = sys.argv[1]

data = np.genfromtxt(fname)
n = data.shape[0]
print("Number of points = " + str(n))

plt.scatter(data[:,0], data[:,1])
plt.grid('on')
plt.show()
