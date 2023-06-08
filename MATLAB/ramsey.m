clear
clc

N = 2;
[Ix,Iy,Iz,IHx,IHy,IHz,sIHz] = prodop([1/2 1/2],[N]);

% Parameters for the NV-centre (in contact with a spin-bath)
hbar = 1;
B0 = 0.055; % [T] Applied DC field
gamma_e = 1.7609e11/2/pi; % Hz Electron gyromagnetic ratio.
gamma_N15 = -4.316e6; % Hz gyromagnetic ratio.
omega_e = gamma_e*B0; % [Hz] Zeeman splitting of electron triplet
omega_N15 = gamma_N15*B0; % [Hz] Zeeman splitting of N14

data = importdata('20221124-1542-40_pulsed_measurement.dat').data;
A_perp = 0; % [Hz] Hyperfine coupling constant, perpendicular to [1, 1, 1] axis. Can be neglected under secular approximation.
A_parallel = 2.16e6; % [Hz] hyperfine coupling constant on the z-axis

n = 500;
rx = expm(-1i*Ix(:,:,1)*pi/2);
ry = expm(-1i*Iy(:,:,1)*pi/2);

rho0 = [1 0]'*[1 0];
rho0 = ry*kron(rho0,eye(2)/2)*ry';

% psize = 100;
% generations = 10;
% inguess = [1e7 1e7 3e7 1];
% poprange = [1e6 0.1e7 1e7 0.8;...
%             5e8 20e7  5e8  3];
% options= optimoptions('ga',...
%     'PopulationSize',psize,...
%     'MaxGenerations',generations,...
%     'PlotFcns',@gaplotbestf,...
%     'Display','final',...
%     'StallGenLimit',generations,...
%     'PopInitRange',poprange,...
%     'InitialPopulation',inguess,...
%     'UseParallel',false,...
%     'StallTimeLimit',inf);    

% [xopt,Q] = ga(@(x) Rabi_opt(x,data,rho0,omega_N15,A_parallel,A_perp,Ix,Iz,Iy,N),4,options);
% off = abs(xopt(1));
% rabi_drive = abs(xopt(2));
% gamma = abs(xopt(3));
% b = abs(xopt(4));
gamma = 4.0267e+07;
b = 0.1;
off = 1.5239e+08;
rabi_drive = 1.1324e+08;
T_init  = data(1,1);
T_final = data(end,1);
t = linspace(T_init,T_final,n);
dt = t(2) - t(1);

H_drive = 2*pi*rabi_drive * Iy(:,:,1);
H_hyp = -2*pi*(off*Iz(:,:,1) + omega_N15*Iz(:,:,2) + (Iz(:,:,1) + 0.5*eye(4))*(A_perp*Ix(:,:,2) + A_parallel*Iz(:,:,2)));
H_total = H_hyp + H_drive;

r = zeros(2^2,2^2,n);
r(:,:,1) = rho0;
r = ry'*r(:,:,1)*ry;
Ez = zeros(1,n);
Ez(1) = trace(r(:,:,1)*Iz(:,:,1));

% dU = @(t) expm(-1i*(H_total)*t);
dU = expm(-1i*(H_total)*dt);

for i = 2:n
    r(:,:,i) = dU*r(:,:,i-1)*dU';
    r(:,:,i) = ry'*r(:,:,i)*ry;
    Ez(i) = real(trace(r(:,:,i)*Iz(:,:,1)))*(exp(-(t(i)*gamma)^b));
end

contrast = abs(data(:,2) - data(:,3)) / (data(:,2) + data(:,3));

figure
% plot(t, Ez, 'b-')
% hold on
plot(t, (Ez-min(Ez))/(max(Ez)-min(Ez)), 'b-')
hold on
plot(data(:,1), (contrast-min(contrast)) / ( max(contrast) - min(contrast) ), 'k-')
% plot(data(:,1), (data(:,2)-min(data(:,2))) / ( max(data(:,2)) - min(data(:,2)) ), 'k-')
% hold on
% plot(data(:,1), (data(:,3)-min(data(:,3))) / ( max(data(:,3)) - min(data(:,3)) ), 'r-')
hold off
xlabel('Time');
ylabel('Prob');

%open quantum transfer