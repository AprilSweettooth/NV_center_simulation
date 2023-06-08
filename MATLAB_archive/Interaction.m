function [omega]=Interaction(r,dir_1, dir_2, s)
% Solve Hamiltonian equation for two NV center coupling.
%
% omega			dipolar coupling (Hz)
% r             distance (nm)
% dir_1         angle for first NV center centered at [r,0,0]
% dir_2         angle for second NV center centered at [0,0,0]
% s_1           spin number for first particle
% s_2           spin number for second particle
%

hbar=1.054571596*10^(-34);								    %reduced Planck's constant (J s)
h=hbar*2*pi;                                                %Planck's constant (J s)
mu_0=1.25663706212*10^(-6);									%magnetic permeapermeability (H m^-1)
mu_B=9.2740100783*10^(-24);                                 %Bohr magneton (J T^-1)
g_e=-2.00231930436256;                                      %electronic g factor

sx = spin_matrix_x(s);                                      %get all spin matrix
sy = spin_matrix_y(s);
sz = spin_matrix_z(s);

dir_1 = dir_1 / norm(dir_1);                                %direction for individual NV center
dir_2 = dir_2 / norm(dir_2);
S1_dot_S2 = sx*sx* (dir_1(1)*dir_2(1)) + sy*sy*(dir_1(2)*dir_2(2)) + sz*sz*(dir_1(3)*dir_2(3));
S1_r = sx*(dir_1(1)*(1));
S2_r = sx*(dir_2(1)*(1));
% disp(S1_r);
% disp(S2_r);
% disp(S1_dot_S2);
Hmatrix = S1_dot_S2 - 3 * S1_r * S2_r;
[coupling]=eig(Hmatrix);                                    %use matlab function eigs to find num_sol eigenfunctions and eigenvalues
r = r*10^(-9);
prefactor = (mu_0 * g_e^2 * mu_B^2) / (4 * pi * h * r^3); 
% disp(coupling);
omega=(prefactor/hbar) * (hbar) * coupling(2);
% disp(omega);
% omega=prefactor;
return