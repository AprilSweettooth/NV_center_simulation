clear
clc
[Ix,Iy,Iz,IHx,IHy,IHz,sIHz] = prodop(1/2,1);
% Parameters for the NV-centre (in contact with a spin-bath)
hbar = 1;
B0 = 0.2359; % [T] Applied DC field
gamma_e = 1.7609e11; % [rad s^-1 T^-1] Electron gyromagnetic ratio.
gamma_N14 = 19.331e6; % [rad s^-1 T^-1 N-14 gyromagnetic ratio.
omega_e = - gamma_e * B0 / (2 * pi); % [Hz] Zeeman splitting of electron triplet
omega_N14 = - gamma_N14 * B0 / (2 * pi); % [Hz] Zeeman splitting of N14
Dzfs = 2.87e9; % [Hz] Zero field splitting of electron triplet.
data=importdata('Rabi1_3.7376GHz_235.9mT.dat').data;

A_perp = 0; % [Hz] Hyperfine coupling constant, perpendicular to [1, 1, 1] axis. Can be neglected under secular approximation.
A_parallel = 2.16e6; % [Hz] hyperfine coupling constant on the z-axis


rabi_drive = 1e6 * 2*pi;
omega = 3.7376e9;
phi = 0;
n = 500;
T_init  = data(1,1);
T_final = data(end,1);
% t = linspace(T_init,T_final,n);
timesteps = 8000;
interval = (T_final-T_init)/timesteps;

sz_e = spin_matrix_z(1);
sx_e = spin_matrix_x(1);
sy_e = spin_matrix_y(1);
Nz_14 = spin_matrix_z(1);
Nx_14= spin_matrix_x(1);
Ny_14 = spin_matrix_y(1);
In = eye(3);
Jz = spin_matrix_z(1);

H_zfs = kron( 2*pi * Dzfs*sz_e*sz_e, In );
H_zee_e = kron( 2*pi * omega_e*sz_e, In );
H_zee_N_14 = kron( In, 2*pi * omega_N14*Nz_14 );
HI_e_N14 = 2*pi* (A_parallel*kron(sz_e,Nz_14));
H_total = H_zfs + H_zee_e;

projection = kron(Jz,In);
gamma = 5*10^8;
b = 1.5;
H_drive = rabi_drive * kron( sx_e, In );

[V0, D] = eig(H_total);
E = diag(D);
[r, c] = find(E==min(E), 1);
vec = V0(:,r);
rho0 = vec * vec';

probs = zeros(1,timesteps);
rho = rho0;

% H = @(t) (H_total + H_drive * exp(-gamma*(t)^b) * cos(2*pi*omega*(t) + phi));
% disp(H_total)
% disp(H(1.3*10^(-9)));
% disp(H(3*10^(-9)));
% disp(data(2,1));
% disp(H(data(10,1))*rho-rho*H(data(10,1)));
for t = 1:timesteps

    H = H_total + H_drive * exp(-(gamma*t*interval)^b) * cos(2*pi*omega*(t*interval) + phi);
    U = expm( - 1i * H * interval * t);
    rho = U*rho0*U';
    % com = -1i*(H(t*interval)*rho - rho*H(t*interval));
    % rho = rho + interval*com;
    % disp(rho);
    p = trace(rho*projection);
    probs(t) = real(p);
    disp(probs(t));
    fprintf("Calculation " + t + " of " + timesteps + "\n");
end


time = T_init:interval:T_final-interval;
plot(time, (probs-min(probs)) / ( max(probs) - min(probs) ))
hold on
plot(data(:,1), (data(:,2)-min(data(:,2))) / ( max(data(:,2)) - min(data(:,2)) ), 'k-')
hold off
xlabel('Time');
ylabel('Prob');