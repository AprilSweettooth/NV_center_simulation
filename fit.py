import matplotlib.pyplot as plot 
import numpy as np

time = np.loadtxt("Rabi1_3.7376GHz_235.9mT.dat",unpack=True)[0]
signal = np.loadtxt("Rabi1_3.7376GHz_235.9mT.dat",unpack=True)[1]
print(signal)
# plot.plot(*np.loadtxt("Rabi1_3.7376GHz_235.9mT.dat",unpack=True), linewidth=2.0)
# plot.show()
