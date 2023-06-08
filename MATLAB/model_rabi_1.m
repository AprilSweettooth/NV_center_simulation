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

data = importdata('Rabi2_1.33GHz_55mT.dat').data;
A_perp = 0; % [Hz] Hyperfine coupling constant, perpendicular to [1, 1, 1] axis. Can be neglected under secular approximation.
A_parallel = 2.16e6; % [Hz] hyperfine coupling constant on the z-axis

n = 500;

rho0 = [1 0]'*[1 0];
rho0 = kron(rho0,eye(2)/2);


psize = 100;
generations = 10;
inguess = [1e7 1e7 3e7 1];
poprange = [1e6 0.1e7 1e7 0.8;...
            5e8 20e7  5e8  3];
options= optimoptions('ga',...
    'PopulationSize',psize,...
    'MaxGenerations',generations,...
    'PlotFcns',@gaplotbestf,...
    'Display','final',...
    'StallGenLimit',generations,...
    'PopInitRange',poprange,...
    'InitialPopulation',inguess,...
    'UseParallel',false,...
    'StallTimeLimit',inf);    

%startpt = [1e6 10e6];
%[x,res] = lsqcurvefit(@Rabi_opt,startpt,xdata,ydata,lb,ub,A,b,Aeq,beq,@nlcon);
[xopt,Q] = ga(@(x) Rabi_opt(x,data,rho0,omega_N15,A_parallel,A_perp,Ix,Iz,Iy,N),4,options);
off = abs(xopt(1));
rabi_drive = abs(xopt(2));
gamma = abs(xopt(3));
b = abs(xopt(4));
% gamma = 4.0267e+08
% b = 3.3165;
% off = 1.5239e+08;
% rabi_drive = 1.1324e+08;
T_init  = data(1,1);
T_final = data(end,1);
t = linspace(T_init,T_final,n);
dt = t(2) - t(1);

H_drive = 2*pi*rabi_drive * Iy(:,:,1);
H_hyp = -2*pi*(off*Iz(:,:,1) + omega_N15*Iz(:,:,2) + (Iz(:,:,1) + 0.5*eye(4))*(A_perp*Ix(:,:,2) + A_parallel*Iz(:,:,2)));
H_total = H_hyp + H_drive;

% projection = kron(Jz,In);


pz = zeros(1,n);
p1 = zeros(1,n);
p0 = zeros(1,n);


dU = expm(-1i*(H_total)*dt);
rho = zeros(2^N,2^N,n);
rho(:,:,1) = rho0;

r1 = TrX(rho(:,:,1),2,[2 2]);
p0(1) = r1(1,1);
p1(1) = r1(2,2);


pz(1) = trace(rho(:,:,1)*Iz(:,:,1));

for i = 2:n
    rho(:,:,i) = dU*rho(:,:,i-1)*dU';
    pz(i) = trace(rho(:,:,i)*Iz(:,:,1));
%     pz(i) = trace(rho(:,:,i)*Iz(:,:,1));
    r1 = TrX(rho(:,:,i),2,[2 2]);
    p0(i) = (real(r1(1,1)))*(exp(-(t(i)/gamma)^b)+0.5);
    p1(i) = r1(2,2);
    fprintf("Calculation " + i + " of " + n + "\n");
end

figure
plot(t, (p0-min(p0))/max((p0)-min(p0)))
hold on
plot(data(:,1), (data(:,2)-min(data(:,2))) / ( max(data(:,2)) - min(data(:,2)) ), 'k-')
% plot(t,p0)
% hold on
% plot(data(:,1), data(:,2))

hold off
xlabel('Time');
ylabel('Prob');

%extract mid value?
%small amp in mid part
%compare 4 var in two data
%ramsey
%open quantum transfer