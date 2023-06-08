import numpy as np
from numpy import pi
from qutip import *
import matplotlib.pyplot as plt
from pylab import *

font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }

sz_e = jmat(1,'z')
sx_e = jmat(1,'x')
sy_e = jmat(1,'y')
Nz_14 = jmat(1,'z')
Nx_14 = jmat(1,'x')
Ny_14 = jmat(1,'y')
Nz_15 = jmat(1/2,'z')
Nx_15 = jmat(1/2,'x')
Ny_15 = jmat(1/2,'y')
In = qeye(3)

# ground_m = basis(3, 0)
# ground_0 = basis(3, 1)
# ground_p = basis(3, 2)

B0 = 0.2359 # [T]
gamma_e = 1.7609e11 # [rad s^-1 T^-1] Electron gyromagnetic ratio.
gamma_N14 = 19.331e16 # [rad s^-1 T^-1 N-14 gyromagnetic ratio.
omega_e = - gamma_e * B0 / (2 * pi) # [Hz] Zeeman splitting of electron triplet
omega_N14 =  gamma_N14 * B0 / (2 * pi) # [Hz] Zeeman splitting of N14
Dzfs = 2.87e9 # [Hz] Zero field splitting of electron triplet.

A_perp = 0 # [Hz] Hyperfine coupling constant, perpendicular to [1, 1, 1] axis. Can be neglected under secular approximation.
A_parallel = 2.16e6 # [Hz] hyperfine coupling constant on the z-axis

rabi_drive = 5e6 * 2*np.pi
omega = 3.7376e9
phi = 0
T           = 10**(-6)
T_init      = 0
T_final     = 50*10**(-6)
timesteps = 10000
# tlist  = np.linspace(0.0, T, timesteps)
tlist  = np.linspace(T_init, T_final, timesteps)

# psi0 = ground_0*ground_0.dag() # Define initial state

H_zfs = tensor( 2*np.pi * Dzfs*sz_e*sz_e, In )
H_zee_e = tensor( 2*np.pi * omega_e*sz_e, In )
# H_zfs = 2*np.pi * Dzfs*sz_e*sz_e
# H_zee_e = 2*np.pi * omega_e*sz_e
H_zee_N_14 = tensor( In, 2*np.pi * omega_N14*Nz_14 )
# H_zee_N_15 = tensor( In, 2*np.pi * omega_e*Nz_15 )
HI_e_N14 = 2*np.pi* (A_parallel*tensor(sz_e,Nz_14))
# H_total = H_zfs + H_zee_e + H_zee_N_14 + HI_e_N14
H_total = H_zfs + H_zee_e

H = Qobj(np.array([ 
    abs(H_total[0][0]), 
    abs(H_total[1][0]),
    abs(H_total[2][0]),
    abs(H_total[3][0]),
    abs(H_total[4][0]),
    abs(H_total[5][0]),
    abs(H_total[6][0]),
    abs(H_total[7][0]),
    abs(H_total[8][0])  ]), dims=[[3,3],[3,3]])
# H = H.reshape([[3,3],[3,3]])
# print(H.eigenstates(sort='low'))

projection = tensor(spin_Jz(1), In)
gamma = 5*10**6
b = 1.5
H_drive = rabi_drive * tensor( sx_e, In )
args = {'w': omega, 'p': phi}
ground_0 = H.eigenstates(sort='low')[1][0]
H = [H, [H_drive, lambda t,args: np.e**(-gamma*t**b)*np.cos(args['w'] * t * 2*np.pi + args['p'])]]
psi0 = ground_0*ground_0.dag()

p = qutip.mesolve(H, projection, tlist, [], [projection], args)
# print(max(np.real(p.expect[0]))) # check prob mag
# print(min(np.real(p.expect[0])))
def NormalizeData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))
normalized_data = NormalizeData(np.real(p.expect[0]))
time = np.loadtxt("Rabi1_3.7376GHz_235.9mT.dat",unpack=True)[0]
signal = np.loadtxt("Rabi1_3.7376GHz_235.9mT.dat",unpack=True)[1]
normalized_signal = NormalizeData(signal)
# H_drive = rabi_drive * sx_e
# args = {'w': omega, 'p': phi}
# ground_0 = H.eigenstates(sort='low')[1][0]
# H = [H, [H_drive, lambda t,args: 2*np.cos(args['w'] * t * 2*np.pi + args['p'])]]
# psi0 = ground_0*ground_0.dag()

p = qutip.mesolve(H, psi0, tlist, [], [psi0], args)
# plt.plot(tlist, normalized_data, 'r', time, normalized_signal, 'b')
plt.plot(tlist, normalized_data)
plt.xlabel("Time", fontdict=font)
plt.ylabel("Occupation probability", fontdict=font)
plt.ylim(0,1)
# plt.legend("ground_0 simulated", "ground_0 experiment")
plt.show()