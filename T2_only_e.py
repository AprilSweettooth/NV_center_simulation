import numpy as np
from numpy import pi
from qutip import *
import matplotlib.pyplot as plt

sz = jmat(3,'z')
sx = jmat(3,'x')
sy = jmat(3,'y')
Ie = qeye(7)

# h = 1.054*(10**(-34))
B0 = 0.007 # [T]
gamma_e = 1.7609*10**11 # [rad s^-1 T^-1] Electron gyromagnetic ratio.
gamma_N14 = 19.331*10**6 # [rad s^-1 T^-1 N-14 gyromagnetic ratio.
omega_e = - gamma_e * B0 / (2 * pi) # [Hz] Zeeman splitting of electron triplet
omega_N14 = - gamma_N14 * B0 / (2 * pi) # [Hz] Zeeman splitting of N14
Dzfs_g = 2.87*10**9 # [Hz] Zero field splitting of electron triplet.
Dzfs_e = 1.42*10**9
Q = - 4.95*10**6 # [Hz] Zero field splitting of N14.

gamma41 = gamma52 = gamma63 = 63e9
gamma47 = 12e9
gamma57 = gamma67 = 80e9
gamma71 = 3.3e9
gamma72 = gamma73 = 2.4e9

H_zfs = tensor( (2*np.pi * (Dzfs_g + Dzfs_e) * sz*sz) )
H_zee_e = tensor( (2*2*np.pi*omega_e*sz) )

H = H_zfs + H_zee_e

A_0 = H.eigenstates(sort='low')[1][0]
A_m = H.eigenstates(sort='low')[1][1]
A_p = H.eigenstates(sort='low')[1][2]
m   = H.eigenstates(sort='low')[1][3]
E_0 = H.eigenstates(sort='low')[1][4]
E_m = H.eigenstates(sort='low')[1][5]
E_p = H.eigenstates(sort='low')[1][6]

c_op = [np.sqrt(gamma41)*(A_0*E_0.dag()), np.sqrt(gamma52)*(A_m*E_m.dag()), np.sqrt(gamma63)*(A_p*E_p.dag()), np.sqrt(gamma47)*(m*E_0.dag()), np.sqrt(gamma57)*(m*E_m.dag()), np.sqrt(gamma71)*(A_0*m.dag()),np.sqrt(gamma72)*(A_m*m.dag()),np.sqrt(gamma73)*(A_p*m.dag())]
# c_op = []
T      = 10**(-9)
timesteps = 1001
tlist  = np.linspace(0.0, T, timesteps)

omega_pulse = 44.2*10**8 * 2*np.pi
phi = 0

d_H = omega_pulse * sx
args = {'w': 1.2e8, 'p': phi}

H = [H, [d_H, lambda t,args: 2*np.cos(args['w'] * t * 2*np.pi + args['p'])]]

result = qutip.mesolve(H, A_0*A_0.dag(), tlist, c_op, [A_0*A_0.dag(), A_m*A_m.dag(), A_p*A_p.dag(), m*m.dag(), E_0*E_0.dag(), E_m*E_m.dag(), E_p*E_p.dag()], args)

# plt.plot(tlist, np.real(result.expect[0]), 'r',  tlist, np.real(result.expect[1]), 'b', tlist, np.real(result.expect[2]), 'y')
plt.plot(tlist, np.real(result.expect[0]), 'r')
plt.xlabel("Time(ns)")
plt.ylabel("Occupation probability")
# plt.legend(("0 groundstate", "-1 groundstate", "+1 groundstate"))
plt.legend("ground state")
plt.show()


# T2 relaxation time for electron triplet system