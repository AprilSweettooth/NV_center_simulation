function [H_tot] = TwoLevelHamiltonian(r, theta, phi)
% excitation from ground state to excited state arXiv:1612.04783v2 [quant-ph] 19 Jun 2017
% plus nucleus interaction & details for excutation parameters
% two level state https://arxiv.org/pdf/1911.04906.pdf
%
%   Constructs the Hamiltonian H of two interacting spins.
%   include the ground state and excited state of each NV
%   

% hbar = 1;
% NV-NV interaction
sigma_1 = zeros(2, 2);
sigma_1(1,1) = cos(theta(1));
sigma_1(1,2) = sin(theta(1))*exp(1i*phi(1));
sigma_1(2,1) = sin(theta(1))*exp(-1i*phi(1));
sigma_1(2,2) = -cos(theta(1));
sigma_2 = zeros(2, 2);
sigma_2(1,1) = cos(theta(2));
sigma_2(1,2) = sin(theta(2))*exp(1i*phi(2));
sigma_2(2,1) = sin(theta(2))*exp(-1i*phi(2));
sigma_2(2,2) = -cos(theta(2));

sigma_prod = zeros(4,4);
sigma_prod(1,1) = 1;
sigma_prod(2,2) = -1;
sigma_prod(3,3) = -1;
sigma_prod(4,4) = 1;
sigma_prod(2,3) = 2;
sigma_prod(3,2) = 2;

mu_0=1.25663706212*10^(-6);									%magnetic permeapermeability (H m^-1)
mu_B=9.2740100783*10^(-24);                                 %Bohr magneton (J T^-1)
g_e=-2.00231930436256;                                      %electronic g factor
prefactor = (mu_0 * g_e^2 * mu_B^2) / (4 * pi); 
H_int = prefactor*(sigma_prod - 3*kron(sigma_1, sigma_2))/r^3;



B0 = 0*513e-4; % [T] Applied DC field
gamma_e = 1.7609e11; % [rad s^-1 T^-1] Electron gyromagnetic ratio.
gamma_N14 = 19.331e16; % [rad s^-1 T^-1 N-14 gyromagnetic ratio.
omega_e = - gamma_e * B0 / (2 * pi); % [Hz] Zeeman splitting of electron triplet
omega_N14 = - gamma_N14 * B0 / (2 * pi); % [Hz] Zeeman splitting of N14
Dzfs = 2870e6; % [Hz] Zero field splitting of electron triplet.
A_perp = 0; % [Hz] Hyperfine coupling constant, perpendicular to [1, 1, 1] axis. Can be neglected under secular approximation.
A_parallel = -2.16e6; % [Hz] hyperfine coupling constant on the z-axis


NV_zfs = [Dzfs, -omega_N14]; % Zero-field splitting of triplet and N respectively.
NV_zeeman = [omega_e, -omega_N14]; % Zeeman interactions
NV_hf_coupling = [A_perp, A_perp, A_parallel];


H = zeros(4); % Initial empty Hamiltonian

% Add S_X operators

sx = spin_matrix_x(0.5);
H = H + NV_hf_coupling(1)*kron(sx', sx);

% Add S_Y operators

sy = spin_matrix_y(0.5);
H = H + NV_hf_coupling(2)*kron(sy', sy);

% Add S_Z operators

sz = spin_matrix_z(0.5);
H = H + NV_hf_coupling(3)*kron(sz', sz);

% Zero-field splitting

for n = 1:length(NV_zfs)
    H = H+ NV_zfs(n)*kron(sz', sz);
end

% Add Zeeman Hamiltonian

for n = 1:length(NV_zeeman)
    H = H+ NV_zeeman(n)*kron(sz', sz);
end

H_tot = 2*H + H_int;
end