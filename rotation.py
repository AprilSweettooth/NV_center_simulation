import numpy as np
from matplotlib import pyplot

import qutip

Omega = 2.87e9 * 2*np.pi
detuning_list  = np.linspace(-2e8, 2e8, 1001)
gamma = np.array([np.sqrt(Omega**2 + detuning**2) for detuning in detuning_list])
T = 10**(-3)

prob = T * 0.5 * Omega**2 * (1-np.sinc(gamma*T)) / gamma**2

# plot the results
pyplot.plot(detuning_list, 1-prob, 'r')
pyplot.xlabel('Time')
pyplot.ylabel('Occupation probability')
pyplot.show()

# FWR approx method to fast the simulation