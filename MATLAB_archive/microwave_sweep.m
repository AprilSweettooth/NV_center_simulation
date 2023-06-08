clear all;
% Parameters for the NV-centre (in contact with a spin-bath)
hbar = 1;
Dzfs = 2870e6; % [Hz] Zero field splitting of electron triplet.
B0 = 0*43e-5; % [T] Applied DC field
gamma_e = 1.7609e11; % [rad s^-1 T^-1] Electron gyromagnetic ratio.
omega_e = - gamma_e * B0 / (2 * pi); % [Hz] Zeeman splitting of electron triplet

phi = pi/4;
f_MW_start = 2.8e9; % [Hz]
f_MW_end = 2.955e9; % [Hz]
calculations = 2;
f_MW_step = (f_MW_end - f_MW_start) / calculations;

spin_e_triplet = 1;

spins = spin_e_triplet; % Spin of the electron triplet and the N-atom (e-e for now)
zfs = [Dzfs, 0]; % Zero-field splitting of triplet and N respectively.
zeeman = [omega_e, 0]; % Zeeman interactions

spin_multiplicities = 2*spins+1;

H_dim = prod(2*spins + 1);
H0 = zeros( H_dim ); % Initial empty Hamiltonian

for n = 1:length(spins)

    % Add S_Z operators
    sz = spin_matrix_z(spins(n));
    % Zero-field splitting
    hzfs = sz*sz*zfs(n);
    % Add Zeeman Hamiltonian
    hzee = sz*zeeman(n);

    H0 = 2*pi * (H0 + hzfs + hzee);
 
end

% Diagonalise the Hamiltonian

[V, D] = eig(H0);

E = diag(D);

[r, c] = find(E==min(E), 1);
vec = V(:,r);

% rho0 = vec * vec';

% Electron triplet basis states

e_basis = { [1;0;0], [0;1;0], [0;0;1] };

% Initialize calculation

avg_probs = zeros(calculations, 1);

delta_H_op = spin_matrix_x(spin_e_triplet);

time_span = 10;
tspan = 0:1:time_span-1;
A = 1;
phase = 0;

rho0 = e_basis{2}*e_basis{2}';

for i = 1:calculations
    
    f_MW = f_MW_start + f_MW_step * i;
    
    probs = zeros(time_span, 1);

    Ht = @(t) A*delta_H_op*sin(f_MW*t + phase);
    
    % [ts, rhos] = solve_tdse(tspan, rho0, @(t) spin_matrix_z(1) + Ht(t));
    rhos = analyse_tdse(tspan, rho0, @(t) H0 + Ht(t));
    for j = 1:time_span
        rho_mat = zeros(3,3);
        k = 1;
        for r = 1:3
            for c = r:3
                rho_mat(r, c) = rhos(j,k);
                rho_mat(c, r) = rhos(j,k)';
                k = k + 1;
            end
        end
        disp(rho_mat);
        probs(j) = vec'*rho_mat*vec;
    end
    % disp(probs);
    avg_probs(i) = real(mean(probs));
    disp(avg_probs(i));
    fprintf("Calculation " + i + " of " + calculations + "\n");
end


% fs = (f_MW_start+f_MW_step):f_MW_step:f_MW_end;
% plot(fs, avg_probs)
% title('Graph of NV-center ESR. Applied DC field: ' + B0 + "T");
% xlabel('Pulse frequency [Hz]');
% ylabel('Prob(T_0)');


function [ts, rhos] = solve_tdse(tspan, rho0, H)

    hbar = 1;
    
    rho_mat = zeros(length(H(0)));
    rho_len = length(rho_mat);

    function drho_dt = odefunc(t, rho)
        drho_dt = zeros(length(rho), 1);
        
        i = 1;
        for r = 1:rho_len
            for c = r:rho_len
                rho_mat(r, c) = rho(i);
                rho_mat(c, r) = rho(i)';
                i = i + 1;
            end
        end
        
        Ht = H(t);
        rho_mat = (- 1i / hbar) * (Ht*rho_mat - rho_mat*Ht);

        i = 1;
        for r = 1:rho_len
            for c = r:rho_len
                drho_dt(i) = rho_mat(r, c);
                i = i + 1;
            end
        end
    end
    
    [ts, rhos] = ode45(@odefunc, tspan, rho0);
    
end

function [rhos] = analyse_tdse(tspan, rho0, H)

    hbar = 1;
    rhos = zeros(1000,length(rho0));
    dt = tspan(2) - tspan(1);
    rho_mat = zeros(length(H(0)));
        
    i = 1;
    for r = 1:3
        for c = r:3
            rho_mat(r, c) = rho0(i);
            rho_mat(c, r) = rho0(i)';
            i = i + 1;
        end
    end
    
    for t = 1:1000
        Ht = H(t);
        rho_mat = rho_mat - dt * (1i / hbar) * (Ht*rho_mat - rho_mat*Ht);
        i = 1;
        for r = 1:3
            for c = r:3
                rhos(t,i) = rho_mat(r, c);
                i = i + 1;
            end
        end
    end
end