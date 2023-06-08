import numpy as np
from numpy import pi
from qutip import *
from qutip.qip.device import Processor
from qutip.bloch import Bloch
import matplotlib.pyplot as plt
from pylab import *
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from qutip.ipynbtools import plot_animation
from qutip.states import qutrit_basis


s11, s22, s33, s12, s23, s31 = qutrit_ops()
is12 = s12 * -1j
is23 = s23 * -1j
sz = s11 - s33
sx = (np.sqrt(2)/2)*(s12 + s12.dag() + s23 + s23.dag())
sy = (np.sqrt(2)/2)*(is12 + is12.dag() + is23 + is23.dag())
I = s11+s22+s33

# print(sx, sy, sz)

minus_one = basis(3, 0)
ground = basis(3, 1)
plus_one = basis(3, 2)

# h = 1.054*(10**(-34))
# B1 = 0.00007 #Max simulable ratio of Bz/By (i.e., B0/B1) is 1000. Time steps in tlist must be increased by factor of 10 accordingly 
B0 = 0.007 # [T]
gamma_e = 1.7609*10**11 # [rad s^-1 T^-1] Electron gyromagnetic ratio.
gamma_N14 = 19.331*10**16 # [rad s^-1 T^-1 N-14 gyromagnetic ratio.
omega_e = - gamma_e * B0 / (2 * pi) # [Hz] Zeeman splitting of electron triplet
omega_N14 = - gamma_N14 * B0 / (2 * pi) # [Hz] Zeeman splitting of N14
Dzfs = 2.87*10**9 # [Hz] Zero field splitting of electron triplet.
P = - 4.95*10**6 # [Hz] Zero field splitting of N14.
omega_pulse = 44.2*10**6 * 2*np.pi
A_perp = 0 # [Hz] Hyperfine coupling constant, perpendicular to [1, 1, 1] axis. Can be neglected under secular approximation.
A_parallel = -2.16e6 # [Hz] hyperfine coupling constant on the z-axis

omega = 2.87e9
phi = 0
T      = 10**(-7)
timesteps = 1001
tlist  = np.linspace(0.0, T, timesteps)

psi0 = ground*ground.dag() # Define initial state
psi1 = plus_one*plus_one.dag() # Define initial state
psi2 = minus_one*minus_one.dag() # Define initial state

H_zfs = Dzfs*sz*sz
# H_zee = omega_e*sz
H_zee = 9.8e6*sz
H0 = 2*np.pi* (H_zfs + H_zee)  # time-independent term
H1 = omega_pulse * sx
args = {'w': omega, 'p': phi}

H = [H0, [H1, lambda t,args: 2*np.cos(args['w'] * t * 2*np.pi + args['p'])]]
# p = qutip.mesolve(H, ground, tlist, [], [psi0, psi1, psi2], args)
p = qutip.mesolve(H, ground, tlist, [], [psi0, psi1, psi2], args)

plt.plot(tlist, np.real(p.expect[0]), 'r',  tlist, np.real(p.expect[1]), 'b', tlist, np.real(p.expect[2]), 'y')
plt.xlabel("Time")
plt.ylabel("Occupation probability")
plt.legend(("0 state", "+1 state", "-1 state"))
plt.show()


# rabi oscillation with microwave drive