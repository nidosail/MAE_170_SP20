%% Given Constants
rho = 8522; % kg/m3
c_s = 380; % J/kgK
k_s = 104; % W/mK
alpha = 3.21e-5; % m2/s
R = .0254; % m
R_err = .0001; % m

%% Calculating T and time values and error
T_s = (mean(volt_al_con')./3).*100;
T_s_err = (std(volt_al_con')./3).*100;
t = mean(time_al_con');
t_err = std(time_al_con');
% various values for calculating theta
T_o = T_s(1);
T_inf = T_s(end);

%% Calculating theta, tau and error
theta_S = (T_s - T_inf)./(T_o - T_inf);
theta_S_err = theta_S.*sqrt((T_s_err./T_s).^2 + (T_s_err(1)./T_o).^2);
tau = (k_s.*t)./(rho.*c_s.*(R.^2));
tau_err = (k_s.*t_err)./(rho.*c_s.*(R.^2));

%% Plotting
y = log(theta_S);
yerr = theta_S_err./theta_S;
x = tau;
xerr = tau_err;
figure(1)
plot(x,y,'r')
hold on
xlabel('\tau  (dimensionless)')
ylabel('ln(\theta_S)')

%% Fitting
[p,S] = polyfit(x(100:500),y(100:500),1);
COVB = (S.R'*S.R)\eye(size(S.R)) *S.normr^2/S.df;
SE = sqrt(diag(COVB));
fprintf('The slope of the linear fit is m = %.4f +/- %.4f \n',p(1),SE(1));
xfit = linspace(x(100),x(500));
yfit = p(1).*xfit + p(2);
plot(xfit,yfit,'b')


