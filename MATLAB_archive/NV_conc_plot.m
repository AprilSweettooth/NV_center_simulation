clear;
close all;

range=40;
npoints=200;
n=logspace(14,22,npoints);

hbar=1.054571596*10^(-34);								    %reduced Planck's constant (J s)
h=hbar*2*pi;                                                %Planck's constant (J s)
mu_0=1.25663706212*10^(-6);									%magnetic permeapermeability (H m^-1)
mu_B=9.2740100783*10^(-24);                                 %Bohr magneton (J T^-1)
g_e=-2.00231930436256;                                      %electronic g factor
prefactor = (mu_0 * g_e^2 * mu_B^2) / (4 * pi * h * (0.55396)^3); 

loglog(n,prefactor*n);			
title('Electronic NV-NV spins');
xlabel('Concentrations of NVs (cm^{-3})') 
ylabel('magnetic dipolar coupling (MHz)') 
