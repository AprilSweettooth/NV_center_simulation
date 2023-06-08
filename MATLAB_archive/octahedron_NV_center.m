function [omega]=octahedron_NV_center(r,center,dir,n)
% Solve Hamiltonian equation for an octahedron shaped NV centers' coupling of total 7 NVs.
%
% omega			dipolar coupling (Hz)
% r             distance (nm)
% The center NV is oriented in vertical y axis and sitting at the origin
%
% center NV orientation (I am not sure if this parameer is redundant as
% rotation is relative
%
% dir consists of 6 individual direction matrix at each corner [ dir_1, dir_2, dir_3, dir_4, dir_5, dir_6 ]
% dir_1         angle for NV center centered at [r,0,0]
% dir_2         angle for NV center centered at [-r,0,0]
% dir_3         angle for NV center centered at [0,r,0]
% dir_4         angle for NV center centered at [0,-r,0]
% dir_5         angle for NV center centered at [0,0,r]
% dir_6         angle for NV center centered at [0,0,-r]
%
% n             number of centers
%

hbar=1.054571596*10^(-34);								    %reduced Planck's constant (J s)
h=hbar*2*pi;                                                %Planck's constant (J s)
mu_0=1.25663706212*10^(-6);									%magnetic permeapermeability (H m^-1)
mu_B=9.2740100783*10^(-24);                                 %Bohr magneton (J T^-1)
g_e=-2.00231930436256;                                      %electronic g factor

sx = spin_matrix_x(1);                                      %get all spin matrix
sy = spin_matrix_y(1);
sz = spin_matrix_z(1);

% connecting_unit = zeros(6,3);
% connecting_unit = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
if n==3
    connecting_unit = [1 0 0; -1 0 0];
elseif n==5
    connecting_unit = [1 0 0; -1 0 0; 0 0 1; 0 0 -1];
end

for i=1:n-1
    dir(i,:) = dir(i,:) / norm(dir(i,:));
end
% dir(1,:) = dir(1,:) / norm(dir(1,:));                                %direction for individual NV center
% dir(2,:) = dir(2,:) / norm(dir(2,:)); 
% dir(3,:) = dir(3,:) / norm(dir(3,:)); 
% dir(4,:) = dir(4,:) / norm(dir(4,:)); 
% dir(5,:) = dir(5,:) / norm(dir(5,:)); 
% dir(6,:) = dir(6,:) / norm(dir(6,:)); 

Hmatrix = zeros(3,3);
for i=1:n-1
    Scenter_dot_Scorner = sx*sx* (center(1)*dir(i,1)) + sy*sy*(center(2)*dir(i,2)) + sz*sz*(center(3)*dir(i,3));
    Scenter_r = sx*(center(1)*connecting_unit(i,1)) + sy*(center(2)*connecting_unit(i,2)) + sz*(center(3)*connecting_unit(i,3));
    Scorner_r = sx*(dir(i,1)*connecting_unit(i,1)) + sy*(dir(i,2)*connecting_unit(i,2)) + sz*(dir(i,3)*connecting_unit(i,3));
    disp(Scenter_dot_Scorner - 3 * Scenter_r * Scorner_r)
    Hmatrix = Hmatrix + (Scenter_dot_Scorner - 3 * Scenter_r * Scorner_r);
%     Hmatix(i) = (Scenter_dot_Scorner - 3 * Scenter_r * Scorner_r);
%     disp(Hmatrix);
%     disp(Scenter_dot_Scorner - 3 * Scenter_r * Scorner_r)
end
% Hmatrix = Hmatrix(1);

[coupling]=eig(Hmatrix);                                    %use matlab function eigs to find num_sol eigenfunctions and eigenvalues
% disp(coupling);
r = r*10^(-9);
prefactor = (mu_0 * g_e^2 * mu_B^2) / (4 * pi * h * r^3); 
% disp(coupling);
omega=(prefactor/hbar) * (hbar) * coupling(2);
% disp(omega);
% omega=prefactor;
return