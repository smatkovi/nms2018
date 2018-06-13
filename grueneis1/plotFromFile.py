from pylab import *
import matplotlib.pyplot as plt

X, Y = [], []
for line in open('EarthOrbit_Euler.dat', 'r'):
  values = [float(s) for s in line.split()]
  X.append(values[0])
  Y.append(values[1])

plt.plot(X, Y)
savefig('figplotFromFile.png')

