% Parameters for the NV-centre (in contact with a spin-bath)
hbar = 1;
B0 = 0*513e-4; % [T] Applied DC field
gamma_e = 1.7609e11; % [rad s^-1 T^-1] Electron gyromagnetic ratio.
gamma_N14 = 19.331e16; % [rad s^-1 T^-1 N-14 gyromagnetic ratio.
omega_e = - gamma_e * B0 / (2 * pi); % [Hz] Zeeman splitting of electron triplet
omega_N14 = - gamma_N14 * B0 / (2 * pi); % [Hz] Zeeman splitting of N14
Dzfs = 2870e6; % [Hz] Zero field splitting of electron triplet.
P = - 4.95e6; % [Hz] Zero field splitting of N14.
% change back to our prefactor


A_perp = 0; % [Hz] Hyperfine coupling constant, perpendicular to [1, 1, 1] axis. Can be neglected under secular approximation.
A_parallel = -2.16e6; % [Hz] hyperfine coupling constant on the z-axis

phi = pi/4;
f_MW_start = 2.8e9; % [Hz]
f_MW_end = 2.955e9; % [Hz]
calculations = 20000;
f_MW_step = (f_MW_end - f_MW_start) / calculations;

tstart = 0;
time_steps = 1000;

spin_e_triplet = 1;
spin_N = 1;

NV_spins = [spin_e_triplet, spin_N]; % Spin of the electron triplet and the N-atom (e-e for now)
% NV_zfs = [Dzfs, -omega_N14]; % Zero-field splitting of triplet and N respectively.
NV_zfs = [0, 0];
% NV_zeeman = [omega_e, -omega_N14]; % Zeeman interactions
NV_zeeman = [0, 0];
% NV_hf_coupling = [A_perp, A_perp, A_parallel];

% Parameters for the spin-bath

% bath_spins = ones(1, 3)*1/2;
% bath_zfs = zeros(1, 3);
% bath_zeeman = zeros(1,3);
bath_spins = [];
bath_zfs = [];
bath_zeeman = [];

% The coupling for the spin-spin interactions of the whole system

num_of_particles = 2 + length(bath_spins);


kx = ones(num_of_particles);
ky = ones(num_of_particles);
kz = ones(num_of_particles);

kx(1,2) = A_perp;
kx(2,1) = A_perp;
ky(1,2) = A_perp;
ky(2,1) = A_perp;
kz(1,2) = A_parallel;
kz(2,1) = A_parallel;

zfs = [NV_zfs bath_zfs]; % Vector containing the zero-field splittings of the whole system
spins = [NV_spins bath_spins]; % Vector containing all the spins of the whole system
zeeman = [NV_zeeman bath_zeeman]; % Vector containing all the zeeman constants for the whole system

spin_multiplicities = 2*spins+1;

H0 = (2*pi)*Tensor_Interaction(kx,ky,kz, spins, zfs, zeeman); % Define unperturbed Hamiltonian

% disp(H0);
% disp(V0);
% disp(E);
[V0, D] = eig(H0);
E = diag(D);
% add in zeeman effect
% add 9 electrons and 9 N14 around single NV center
% add spin-spin interaction coefficient

% (bath spins) - random magnetic field
% master equations for time evolution with density matrix
% time evolution time dependence rabi pulse/coherent drive


[r, c] = find(E==min(E), 1);
vec = V0(:,r);
% what is find func?
rho0 = vec * vec';

% Electron triplet basis states

e_basis = { [1;0;0], [0;1;0], [0;0;1] };

% Initialize calculation

omega_1 = 7e6;

zfs_e_op = kron(spin_matrix_z(1) * spin_matrix_z(1), kron_id_chain(spin_multiplicities(2:length(spin_multiplicities))));
sz_op = kron(spin_matrix_z(1), kron_id_chain(spin_multiplicities(2:length(spin_multiplicities))));

Hx = (1 / sqrt(2))* [ 0 1 0 ; 1 0 1 ; 0 1 0];
Hx = kron(Hx, kron_id_chain(spin_multiplicities(2:length(spin_multiplicities))));

Hy = (1 / sqrt(2))* [ 0 -1i 0 ; 1i 0 1i ; 0 -1i 0];
Hy = kron(Hy, kron_id_chain(spin_multiplicities(2:length(spin_multiplicities))));


avg_probs = zeros(calculations, 1);

for i = 1:calculations
    
    f_MW = f_MW_start + f_MW_step * i;
    %f_MW = 1.431e9;
    % Add perturbation

    H_e = - f_MW * zfs_e_op;

    H = H0 + (2*pi)*( H_e + (omega_1 / sqrt(2))*(cos(phi) * Hx + sin(phi) * Hy) );
    
    rho = rho0;
    probs = zeros(time_steps, 1);
    
    tend = 2e3/f_MW;
    delta_t = (tend - tstart) / time_steps;
    
    for t = 0:(time_steps-1)
        U = expm( - 1i * H * delta_t * t);
        rho = U*rho0*U';
        
        % Partial trace
        
        rho_A = zeros(3);
        
        for j = 1:3
            for k = 1:3
                for l=1:3
                    rho_A = rho_A + rho( 3*(j-1) + l , 3*(k - 1) + l ) * ( e_basis{j} * e_basis{k}' );
                end
            end
        end
        
        p = e_basis{2}' * rho_A * e_basis{2};
        
        probs(t+1) = p;
    end
    

    avg_probs(i) = real(mean(probs));
    disp(avg_probs(i));
    fprintf("Calculation " + i + " of " + calculations + "\n");
end


fs = (f_MW_start+f_MW_step):f_MW_step:f_MW_end;
plot(fs, avg_probs)
% save_dat = [fs' avg_probs];
% filename = "./results/matlab_nv_sim" + datestr(now,'mm.dd.yyyy-HH.MM.SS');
% save(filename + ".dat", "save_dat");
% title('Graph of NV-center ESR. Applied DC field: ' + B0 + "T");
xlabel('Pulse frequency [Hz]');
ylabel('Prob(T_0)');
% pl = plot(fs, avg_probs);
% saveas(pl,char(filename + ".pdf"));
% saveas(pl,char(filename + ".fig"));



function res = kron_id_chain( id_dims )
    %KRON_ID_CHAIN Creates a chain of tensor products of identity matrices of
    %specified dimension.

    res = eye(id_dims(1));

    for k = 2:length(id_dims)
        res = kron(res, eye(id_dims(k)));
    end
end