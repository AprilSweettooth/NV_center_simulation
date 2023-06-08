%-----------------------------------------------------------------------%
%------------- Rabi flopping for a two level atom-----------------------%
%-----------------------------------------------------------------------%

N=1000; % no of time grid points
t=linspace(0,0.25e-6,N); % time array in micro seconds
dt=t(2)-t(1); % step size in time
%---- initial conditions-------------
cg(1:N)=1; % initial ground state 
ce(1:N)=0; % initial excited state
Omega_rabi=2*pi*10*1e6; % Rabi frequency 2pi*10 Mhz
delta=[0 Omega_rabi 2*Omega_rabi]; % different detuning frequencies
for j=1:3 % loop for different detuning frequencies
    
    
for i=1:N-1 % finite difference loop for advancing in time
    
    cg(i+1)=-(1i*ce(i)*Omega_rabi*dt*(exp(1i*delta(j)*t(i))/2))+cg(i);
    ce(i+1)=-(1i*cg(i)*Omega_rabi*dt*(exp(-1i*delta(j)*t(i))/2))+ce(i);
    
end
Cg(j,:)=cg;
Ce(j,:)=ce;
end
% plotting the probabilities in ground states spins
figure(1)
subplot(3,1,1)
plot(t,abs(Cg(1,:).*Cg(1,:)),'linewidth',2)
hold on
plot(t,abs(Ce(1,:).*Ce(1,:)),'r','linewidth',2)
axis([0 t(N) 0 1.1])
h=gca; 
get(h,'FontSize') 
set(h,'FontSize',14)
axis([0 t(N) 0 1.1])
fh = figure(1);
set(fh, 'color', 'white'); 
subplot(3,1,2)
plot(t,abs(Cg(2,:).*Cg(2,:)),'linewidth',2)
hold on
plot(t,abs(Ce(2,:).*Ce(2,:)),'r','linewidth',2)
ylabel('|cg^2| and |ce^2| probabilities','fontSize',14);
axis([0 t(N) 0 1.1])
subplot(3,1,3)
plot(t,abs(Cg(3,:).*Cg(3,:)),'linewidth',2)
hold on
plot(t,abs(Ce(3,:).*Ce(3,:)),'r','linewidth',2)
hold on
axis([0 t(N) 0 1.1])
xlabel('time','fontSize',14);