ground_0 = [1 0 0]'; % Quantum ground state 0 = [1 0 0]
ground_1 = [0 1 0]'; % Quantum ground state 1 = [0 1 0]
excited = [0 0 1]'; % Quantum excited state = [0 0 1]

dim = 3; % Dimension of the total Hilbert space
Is = eye ( dim ) ; % Identity matrix for the total Hilbert space
hbar = 1; % reduced planck constant


initial_rho = ground_0 * ground_0'; % Initial density matrix projection
Nt = 3000; % Number of steps for time
ti = 0; % Intial time
tf = 30 ; % Final time
dt = (tf - ti ) /( Nt -1) ; % Step for time
t = ti : dt : tf ; % Time vector
c = zeros(size(t));
rho = initial_rho; % Initial density matrix

% for n =1: length (t) % Iteration to find general sulution of rho (t)
%     T = n*dt;
%     H = One_NV_Hamiltonian(T);
%     rho = rho - 1i*dt*(H*rho - rho*H); % General sulution for rho (t)
%     c(n) = ground_0'*rho*ground_0;
% end

for n =1: length (t) % Iteration to find general sulution of rho (t)
    T = n*dt;
    H = One_NV_Hamiltonian(T);
    rho = rho - 1i*dt*(H*rho - rho*H); % General sulution for rho (t)
    c(n) = ground_0'*rho*ground_0;
end

figure ()
hold on
plot (t , real (c))
hold off
xlabel ('time')
ylabel ('Probability')

function H = One_NV_Hamiltonian(t)  % Hamiltonian of the system

ground_0 = [1 0 0]'; % Quantum ground state 0 = [1 0 0]
ground_1 = [0 1 0]'; % Quantum ground state 1 = [0 1 0]
excited = [0 0 1]'; % Quantum excited state = [0 0 1]

omega_l = 16e6; % Rabi frequency l
G0 = 1; % energy of ground state 0
G1 = 1; % energy of ground state 1
E1 = 1; % energy of ground excited state
E = 1; % laser electric field
d13 = 1; % transition dipole moment
Omega = E*d13; % Rabi frequency of layser electric field at omega_l
lambda = 1; % Rabi frequency of ESR transition between m0 and m1 at omega_m
omega_m = 1; % Rabi frequency m
H = G0*kron(ground_0, ground_0') + G1*kron(ground_1, ground_1') + E1*kron(excited, excited') - Omega*cos(omega_l*t)*(kron(ground_0, excited')+kron(excited, ground_0')) - lambda*cos(omega_m*t)*(kron(ground_0, ground_1')+kron(ground_1, ground_0'));
end

function H = Two_NV_Hamiltonian(t)

ground_0 = [1 0 0]'; % Quantum ground state 0 = [1 0 0]
ground_1 = [0 1 0]'; % Quantum ground state 1 = [0 1 0]
excited = [0 0 1]'; % Quantum excited state = [0 0 1]

omega_l = 16e6; % Rabi frequency l
G0 = 1; % energy of ground state 0
G1 = 1; % energy of ground state 1
E1 = 1; % energy of ground excited state
E = 1; % laser electric field
d13 = 1; % transition dipole moment
Omega = E*d13; % Rabi frequency of layser electric field at omega_l
lambda = 1; % Rabi frequency of ESR transition between m0 and m1 at omega_m
omega_m = 1; % Rabi frequency m
H = G0*kron(ground_0, ground_0') + G1*kron(ground_1, ground_1') + E1*kron(excited, excited') - Omega*cos(omega_l*t)*(kron(ground_0, excited')+kron(excited, ground_0')) - lambda*cos(omega_m*t)*(kron(ground_0, ground_1')+kron(ground_1, ground_0'));
end