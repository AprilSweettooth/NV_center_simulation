function [Ix,Iy,Iz,IHx,IHy,IHz,sIHx,sIHy,sIHz] = prodop1(spinNumbers,spinList)
% *****************************
% NMRQC-GROUP, IISER Pune, 2015
% *****************************
if nargin < 2; spinList = ones(length(spinNumbers),1); end;
M = length(spinList);
N = sum(spinList);
spins = [];
for k = 1:M
  spins = [spins spinNumbers(k)*ones(1,spinList(k))];
end
D = prod(2*spins+1);
Ix = zeros(D,D,N); Iy = Ix; Iz = Ix; 

id = eye(N,N);
for k=1:N
  Px=1; Py=1; Pz=1;
  for j=1:N
    [pe,px,py,pz] = genBasicOp(spins(j));
    Px = kron(Px,(id(k,j)*px + (1-id(k,j))*pe));
    Py = kron(Py,(id(k,j)*py + (1-id(k,j))*pe));
    Pz = kron(Pz,(id(k,j)*pz + (1-id(k,j))*pe));
  end
  Ix(:,:,k) = Px;
  Iy(:,:,k) = Py;
  Iz(:,:,k) = Pz;
end
firstsp = 1;
for k = 1:M
  lastsp = firstsp + spinList(k) - 1;
  IHx(:,:,k) = sum(Ix(:,:,firstsp:lastsp),3);
  IHy(:,:,k) = sum(Iy(:,:,firstsp:lastsp),3);
  IHz(:,:,k) = sum(Iz(:,:,firstsp:lastsp),3);
  firstsp = lastsp + 1;
end
sIHx = sum(IHx,3);
sIHy = sum(IHy,3);
sIHz = sum(IHz,3);
function [Pe,Px,Py,Pz] = genBasicOp(j)
Px = 0; Py = 0; Pz = 0;
m = j:-1:-j;   
Pz = diag(m);  Pe = eye(size(Pz));
Pp = zeros(size(Pz));  Pm = Pp;  
for k = 2:length(m); Pp(k-1,k) = sqrt((j-m(k))*(j+m(k)+1)); end;
for k = 1:length(m)-1; Pm(k+1,k) = sqrt((j+m(k))*(j-m(k)+1)); end;
Px = (Pp + Pm)/2;  Py = (Pp - Pm)/2i;

