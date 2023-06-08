clear all
% from igor Iovchinsky phd thesis
mu_B=9.2740100783*10^(-24); %Bohr magneton (J T^-1)
g_e=-2.00231930436256; %electronic g factor
prefactor = -0.5*1i * mu_B * g_e; %prefactor for function S


N=3000; % no of time grid points
t=linspace(0,0.25e-7,N); % time array in micro seconds
dt=t(2)-t(1); % step size in time
beta = zeros(size(t));


% for n =1: N % Iteration to find general sulution of rho (t)
%     T = n*dt;
%     fun = @(x,lambda, omega) prefactor*exp(1i*(lambda+1*omega)*x);
%     c = integral(@(x) fun(x,2.87*10^9,10^9),0,T);
%     beta(n) = abs(c)^2;
% end


%---- initial conditions-------------
cg(1:N)=1; % initial ground state 
ce(1:N)=0; % initial excited state
Omega=10^9; % Rabi frequency 1 Ghz
delta = 2.87*10^9;
omega=[0 delta 2*delta]; % different detuning frequencies
for j=1:3 % loop for different detuning frequencies
    upsilon = omega(j) - delta;
    gamma = sqrt(upsilon^2 + Omega^2);
for i=1:N-1 % finite difference loop for advancing in time
    cg(i) = ((Omega/gamma) * sin(gamma*dt*i*0.5))^2;
    ce(i) = 1-cg(i);
    
end
Cg(j,:)=cg;
Ce(j,:)=ce;
end
% plotting the probabilities in ground states spins
figure(1)

subplot(3,1,1)
plot(t,Cg(1,:),'linewidth',2)
hold on
plot(t,Ce(1,:),'linewidth',2)
axis([0 t(N) 0 1.1])
h=gca; 
get(h,'FontSize') 
set(h,'FontSize',14)
axis([0 t(N) 0 1.1])
fh = figure(1);
set(fh, 'color', 'white'); 

subplot(3,1,2)
plot(t,Cg(2,:),'linewidth',2)
hold on
plot(t,Ce(2,:),'linewidth',2)
ylabel('|cg^2| and |ce^2| probabilities','fontSize',14);
axis([0 t(N) 0 1.1])

subplot(3,1,3)
plot(t,Cg(3,:),'linewidth',2)
hold on
plot(t,Ce(3,:),'linewidth',2)
axis([0 t(N) 0 1.1])
xlabel('time','fontSize',14);