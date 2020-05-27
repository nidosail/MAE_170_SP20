clear all
close all
clc
% [C2H4 O2 C2H4O CO2 H2O] = [A B C D E] = [1 2 3 4 5]
rho_cat = 890; % kgcat/m3
phi = .5;
rho_bed = (1-phi)*rho_cat;
ntube = 5523;
l = 8.7; % m
nT0 = .89*ntube;
n0 = [.25 .08 0 .07 0 .47 .13].*nT0;
T0 = 463; % K
P0 = 21.7; % bar
Y0 = [n0 T0 P0];
%% Single D
% d = 5;
% r = d./2;
% V = pi*r.^2*l*ntube;
% Wt = rho_bed*V;
% Wspan = [0 Wt];
% [W, Y] = ode15s(@(W,Y) fun(W,Y,d),Wspan,Y0);
% figure(1)
% plot(W,Y(:,9)./P0)
% figure(2)
% plot(W,Y(:,1).*44.*3600.*24.*1e-6)
%% First plots
d = linspace(1, 20);
subplot(2,2,1)
for i = 1:100
    D = d(i);
    r = D/2;
    V = pi*r.^2*l*ntube;
    Wt = rho_bed*V;
    Wspan = [0 Wt];
    [W,Y] = ode45(@(W,Y) fun(W,Y,D),Wspan,Y0);
    EOvec(i) = Y(end,3);
    Pvec(i) = Y(end,9);
    Tvec(i) = Y(end,8);

end

subplot(2,2,1)
plot(d,EOvec.*44.*3600.*24.*1e-6)
% plot(W,Y(:,1))
xlabel('Diameter of tube')
ylabel('Flow rate of C_2H_4O (mton/day)')

subplot(2,2,2)
plot(d,Pvec./P0)
title('pressure drop')
subplot(2,2,3)
plot(d,Tvec)

% function for ode45
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
ntot = [nA nB nC nD nE nF nG]; % mol
y = [nA nB nC nD nE nF nG]./sum(ntot); % AU
pP = y.*P; % percent
%  Constants
R = 8.314; 
Cp = [63 31 75 45 35 47 30]; % J/molK
dHf = [52500 0 -52600 -393500 -241800 -74500 0]; % J
mm = [28 32 44 44 18 16 28]; % g/mol
ntube = 5523;
nT0 = .89*ntube; % mol/s
n0 = [.25 .08 0 .07 0 .47 .13].*nT0; % mol/s
%   Catalyst and reactor
rho_cat = 890; % kgcat/m3
phi = .5; % void fraction
rho_bed = (1-phi)*rho_cat; % kgbed/m3
d_p = 7.74e-3; % m
%   kinetic parameters
A = [1.33e5 1.8e6]; % mol/kgs
E_a = [60.7 73.2].*1000; % J
K_E = [6.5 4.33]; % 1/bar
n = [.58 .3];
%  enthalpy
dCp = [nu(1,:)*Cp' nu(2,:)*Cp']; % J/molK
dH1 = nu(1,:)*dHf' + dCp(1)*(T - 298); % J/mol
dH2 = nu(2,:)*dHf' + dCp(2)*(T - 298);
dH = [dH1 dH2];
%  rates
k = A.*exp(-E_a./(R.*T)); % mol/kgcat
r = (k.*pP(1).*pP(2).^n)./((1 + K_E.*pP(1)).^2); % mol/kgcat-s

rp = r./rho_bed;
rA = r*nu(:,1);
rB = r*nu(:,2);
rC = r*nu(:,3);
rD = r*nu(:,4);
rE = r*nu(:,5);
rF = 0;
rG = 0;
%  Re and fp
mT = (sum(mm.*n0))/1000; % kg/s
mu = 10^-5;
Re = (4*mT)/(pi*d*mu); % unitless
fp = (160 + 3*(Re/(1-phi))^0.9)*((1-phi)^2/(phi^3*Re));
v = (sum(ntot)*8.314e-5*T)/P; % m3
%  Heat Transfer
rU = (1/60) + (1/1500);
U = 1/rU;
a = 4/d; % 1/m
T_u = 503; % K
dQp = (U*a*(T_u - T))/rho_bed;
%  Material Balance
dnAdW = rA/rho_bed;
dnBdW = rB/rho_bed;
dnCdW = rC/rho_bed;
dnDdW = rD/rho_bed;
dnEdW = rE/rho_bed;
dnFdW = rF/rho_bed;
dnGdW = rG/rho_bed;
%  Pressure drop
dPdW = (-v*fp*64*mT)/(pi^3*d^6*d_p*rho_bed);
%  Energy Balance
dTdW = (dQp - sum(rp.*dH))/(sum(ntot.*Cp));
dYdW = [dnAdW; dnBdW; dnCdW; dnDdW; dnEdW; dnFdW; dnGdW; dTdW; dPdW];
end