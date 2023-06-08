% define two NV parameters
r = 10^-9;
theta = [0 0];
phi = [pi 0];
H = TwoLevelHamiltonian(r, theta, phi)*10^-9; % Hamiltonian of the system
% initial density matrix
[V, D] = eig(H);
E = diag(D);
r = find(E==min(E), 1);
vec = V(:,r);
rho = vec * vec';
initial_rho = rho;

Nt = 3000; % Number of steps for time
ti = 0; % Intial time
tf = 4; % Final time
dt = (tf - ti ) /( Nt -1) ; % Step for time
t = ti : dt : tf ; % Time vector
c = zeros(size(t));


for n =1: length (t) % Iteration to find general sulution of rho (t)
    T = n*dt;
    H = TwoLevelHamiltonian(r, theta, phi);
    rho = rho - 1i*dt*(H*rho - rho*H); % General solution for rho (t)
    c(n) = trace(initial_rho*rho);
end
% probs = zeros(Nt, 1);
% e_basis = { [1;0;0;0], [0;1;0;0], [0;0;1;0] [0;0;0;1]};
% for t = 0:(Nt-1)
%     U = expm( - 1i * H * dt * t);
%     rho = U*rho0*U';
%     rho_A = zeros(4);
%     
%     for j = 1:4
%         for k = 1:4
%             for l=1:4
%                 rho_A = rho_A + rho( 4*(j-1) + l , 4*(k - 1) + l ) * ( e_basis{j} * e_basis{k}' );
%             end
%         end
%     end
%     
%     p = e_basis{2}' * rho_A * e_basis{2};
%     
%     probs(t+1) = p;
% end

figure ()
hold on
plot (t , real (c))
hold off
xlabel ('time')
ylabel ('Probability')