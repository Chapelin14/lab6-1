#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

x1 = np.loadtxt("x1.dat")
x2 = np.loadtxt("x2.dat")
data = np.loadtxt("direction.dat")
Z = np.loadtxt("f.dat")
X1, X2 = np.meshgrid(x1, x2)
plt.contour(X1, X2, Z)
plt.plot(data[:,0], data[:,1])
plt.show()
