function [Q] = Rabi_opt(x,data,rho0,omega_N15,A_parallel,A_perp,Ix,Iz,Iy,N)

off = abs(x(1));
rabi_drive = abs(x(2));
gamma = abs(x(3));
b = abs(x(4));
t = data(:,1);
dt = t(2) - t(1);
n = length(t);

H_drive = 2*pi*rabi_drive * Iy(:,:,1);
H_hyp = -2*pi*(off*Iz(:,:,1) + omega_N15*Iz(:,:,2) + (Iz(:,:,1) + 0.5*eye(4))*(A_perp*Ix(:,:,2) + A_parallel*Iz(:,:,2)));
H_total = H_hyp + H_drive;


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
    p0(i) = real(r1(1,1))*(exp(-(t(i)/gamma)^b)+0.5);
    p1(i) = r1(2,2);
end

p0 = (p0-min(p0))/max((p0)-min(p0));
d = (data(:,2)-min(data(:,2))) / (max(data(:,2)) - min(data(:,2)));
Cof = (p0-d').^2;
 
Q = sum(Cof);
Q = real(Q);
disp(Q);
if Q > 0
    return
end
end
