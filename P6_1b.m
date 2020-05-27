clear all
close all
clc
% [C2H4 O2 C2H4O CO2 H2O] = [A B C D E] = [1 2 3 4 5]
rho_cat = 890; % kgcat/m3
phi = .5;
rho_bed = (1-phi)*rho_cat;
ntube = 5523;
l = 8.7; % m
d = 16.36; % m
r = d./2;
V = pi*r.^2*l*ntube;
W = rho_bed*V;
nT0 = .89*ntube;
n0 = [.25 .08 0 .07 0 .47 .13].*nT0;
T0 = 463; % K
P0 = 21.7; % bar
Wspan = [0 W(end)];
Y0 = [n0 T0 P0];
D = d;
[W,Y] = ode45(@(W,Y) fun(W,Y,D),Wspan,Y0);

subplot(2,2,1)
plot(W,Y(:,3))
hold on
plot(W,Y(:,1))
plot(W,Y(:,2))
plot(W,Y(:,4))
plot(W,Y(:,5))
hold off
xlabel('Mass of catalyst')
ylabel('Molar flow rates)')

subplot(2,2,2)
Pdrop = Y(:,9)./P0;
plot(W,Pdrop)
title('pressure drop')
subplot(2,2,3)
plot(W,Y(:,8))

function dYdW = fun(W,Y,d)
nA = Y(1);
nB = Y(2);
nC = Y(3);
nD = Y(4);
nE = Y(5);
nF = Y(6);
nG = Y(7);
T = Y(8);
P = Y(9);
nu = [-1 -1/2 1 0 0 0 0; -1 -3 0 2 2 0 0];
ntot = [nA nB nC nD nE nF nG];
y = [nA nB nC nD nE nF nG]./sum(ntot);
pP = y.*P;
%  Constants
R = 8.314; 
Cp = [63 31 75 45 35 47 30];
dHn = [52500 0 -52600 -393500 -241800 -74500 0];
mm = [28 32 44 44 18 16 28];
ntube = 5523;
nT0 = .89*ntube;
n0 = [.25 .08 0 .07 0 .47 .13].*nT0;
%   Catalyst and reactor
rho_cat = 890; % kgcat/m3
phi = .5;
rho_bed = (1-phi)*rho_cat;
d_p = 7.74e-3; % m
%   kinetic parameters
A = [1.33e5 1.8e6];
E_a = [60.7 73.2].*1000;
K_E = [6.5 4.33];
n = [.58 .3];
%  enthalpy
dH1 = nu(1,:)*dHn';
dH2 = nu(2,:)*dHn';
dH = [dH1 dH2];
%  rates
k = A.*exp(-E_a./(R.*T));
r = (k.*pP(1).*pP(2).^n)./((1 + K_E.*pP(1)).^2);
rp = r./rho_bed;
rA = r*nu(:,1);
rB = r*nu(:,2);
rC = r*nu(:,3);
rD = r*nu(:,4);
rE = r*nu(:,5);
rF = 0;
rG = 0;
%  Re and fp
mT = (mm*n0')/1000;
mu = 10^-5;
Re = 4*mT/(pi*d*mu);
fp = (150 + 1.75*(Re/(1-phi))*((1-phi)^2/(phi^3*Re)));
v = sum(ntot)*8.314e-5*T/P;
%  Heat Transfer
dQ = 0;
%  Material Balance
dnAdW = rA/rho_bed;
dnBdW = rB/rho_bed;
dnCdW = rC/rho_bed;
dnDdW = rD/rho_bed;
dnEdW = rE/rho_bed;
dnFdW = rF/rho_bed;
dnGdW = rG/rho_bed;
%  Pressure drop
dPdW = -v*fp*64*mT/(pi^3*d^6*d_p);
%  Energy Balance
dTdW = (dQ - rp*dH')/(ntot*Cp');
dYdW = [dnAdW; dnBdW; dnCdW; dnDdW; dnEdW; dnFdW; dnGdW; dTdW; dPdW];
end