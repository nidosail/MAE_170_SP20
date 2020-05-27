clear all;
close all;
clc
Tguess = 434;
T0 = Tguess + 273;
Tu0 = Tguess + 273;
P0 = 10;
v0 = 20; % L/s
R = .08206;
nA0 = (P0*v0)/(R*T0);
n0 = [nA0 0 0];
Y0 = [n0 T0 P0 Tu0];
Wrange = [0 80];
[W, Y] = ode45(@fun, Wrange, Y0);
fprintf('Tu(W) = %.2f C \n',Y(end,end)-273);
disp(Y(end,4)-273)
function dYdW = fun(W,Y)
nA = Y(1);
nB = Y(2);
nC = Y(3);
T = Y(4);
P = Y(5);
Tu = Y(6);
n = [nA nB nC];
%  Constants
d = .0508; % m
E_a = 31400; % J/mol
k_450 = .133; % L/kgs
rhobed = 900/1000; % kg/L
Cp = [40 25 15]; % J/molC
dHf = [-70000 -50000 -40000]; % J/mol
v0 = 20; % L/s
P0 = 10;
RaL = .08206;
v = (sum(n)*RaL*T)/P0;
nu = [-1 1 1];
R = 8.3145;
Tu0 = 383;
%  Calculations
Ca = nA/v;
k = k_450*exp((E_a/R)*((1/450) - (1/T)));
r = k*Ca;
dCp = nu*Cp';
dH = nu*dHf' + dCp*(T-298);
nA0 = (P0*v0)/(RaL*Tu0);
%  rates
rp = r/rhobed;
rA = nu(1)*r;
rB = nu(2)*r;
rC = nu(3)*r;
%  Heat transfer calcs
Cpu = Cp(1);
nU = nA0;
U = 30;
a = 4/d;
dQ = (U*a/rhobed)*(Tu - T);
%  Material Balance
dnAdW = rA/rhobed;
dnBdW = rB/rhobed;
dnCdW = rC/rhobed;
%  Energy Balance
dTdW = (dQ - rp*dH)/(n*Cp');
dTudW = dQ/(nU*Cpu);
dPdW = 0;
dYdW = [dnAdW; dnBdW; dnCdW; dTdW; dPdW; dTudW];
end
