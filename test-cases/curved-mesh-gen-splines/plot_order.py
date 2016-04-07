#! /usr/bin/env python3

""" Plots an order-of-accuracy plot using data points given by a command line argument.
 Uses linear regression to compute slope of the data.
"""

import sys
import numpy as np
from matplotlib import pyplot as plt

fname = sys.argv[1]

data = np.genfromtxt(fname)
n = data.shape[0]
print("Number of points = " + str(n))
data = np.log10(data)

psigy = data[:,1].sum()
sigx = data[:,0].sum()
sigx2 = (data[:,0]*data[:,0]).sum()
psigxy = (data[:,1]*data[:,0]).sum()

pslope = (n*psigxy-sigx*psigy)/(n*sigx2-sigx**2)

plt.plot(data[:,0],data[:,1],'o-')
plt.title("Grid-refinement Study: slope = "+str(pslope))
plt.xlabel("Log mesh size")
plt.ylabel("Log error")
plt.show()
