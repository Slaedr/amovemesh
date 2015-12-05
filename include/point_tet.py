import numpy as np
from numpy.linalg import det

v = np.array([[-2,3,-2],[6,3,-2],[-2,3,6],[-2,-5,-2]],dtype=np.float64)
pt = np.array([0,0,0])

eqnmat = np.zeros((5,5),dtype=np.float64)

eqnmat[0][0] = pt[0]**2 + pt[1]**2 + pt[2]**2
eqnmat[0][1] = pt[0]
eqnmat[0][2] = pt[1]
eqnmat[0][3] = pt[2]
eqnmat[0][4] = 1

for i in range(0,4):
	eqnmat[i+1][0] = v[i][0]**2 + v[i][1]**2 + v[i][2]**2
	eqnmat[i+1][1] = v[i][0]
	eqnmat[i+1][2] = v[i][1]
	eqnmat[i+1][3] = v[i][2]
	eqnmat[i+1][4] = 1.0

val = det(eqnmat)
print("Value is " + str(val))
