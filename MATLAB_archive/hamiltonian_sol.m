range=40;
npoints=200;
r=linspace(1,range,npoints);

% For two NV centers
dir_1 = [1,1,0];
dir_2 = [-1,0,0];
dir3 = [0 1 0; 0 1 0];
dir4 = [0 1 0; 0 1 0; 0 1 0; 0 1 0];
for i=1:npoints
     [omega2(i)]=two_nv_Hamiltonian(r(i),dir_1, dir_2);
     [omega3(i)]=octahedron_NV_center(r(i), [0,1,0], dir3, 3);
     [omega4(i)]=octahedron_NV_center(r(i), [0,1,0], dir4, 5);
end
x=linspace(1,range+1,npoints);
semilogy(x,omega2);			
hold on
semilogy(x,omega3);	
hold on
semilogy(x,omega4);	
legend('2 centers','3 centers','4 centers');
title('NV center coherence time comparison');
xlabel('Distance r (nm)') 
ylabel('Dipolar coupling \Omega_{dip} (Hz)') 
% saveas(gcf,'multicenter_comparison.png')
% dir = [0 1 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0];
% for i=1:npoints
%     [omega(i)]=octahedron_NV_center(r(i), [0,1,0], dir);
% end
% x=linspace(1,range+1,npoints);
% semilogy(x,omega);	
% legend('4 centers');
% title('octahedron NV center coherence time')
% xlabel('Distance r (nm)') 
% ylabel('Dipolar coupling \Omega_{dip} (Hz)') 