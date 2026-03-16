% file: simu_benchmark_se.m
%% Three-State DMP with Time-Limited UI and Endogenous Search Effort
% Unknowns: (theta, sI, sU, w)
% Equilibrium system 

clearvars; close all; clc;

cd(fileparts(mfilename('fullpath')))
addpath('Toolbox')

kpic = 1;

%% SECTION 1: Parameters
r     = 0.05;
beta  = 0.5;
delta = 0.15;

kappa = 0.15;   % vacancy posting cost 
y     = 1;

% Matching: f(theta)=A*theta^(1-mu), q(theta)=A*theta^(-mu)
mu = 0.5;      
A  = 0.75;     

% Policy
bI     = 0.25;
bU     = 0.10;
lambda = 0.50;

% Effort cost: c(s)=0.5*c0*s^2
c0 = 1.0;

% Pack parameters
par.r = r; par.beta = beta; par.delta = delta;
par.kappa = kappa; par.y = y;
par.mu = mu; par.A = A;
par.bI = bI; par.bU = bU; par.lambda = lambda;
par.c0 = c0;

%% SECTION 2: Solve for (theta, sI, sU, w)
options = optimset('display','iter', ...
                   'Algorithm','trust-region-dogleg', ...
                   'MaxIter',2000, ...
                   'Tolfun',1e-10, ...
                   'TolX',1e-10);

% Initial guess: [theta, sI, sU, w]
x0 = [0.5, 0.5, 0.5, 0.5];

[x, ~, exitflag] = fsolve(@(x) equilibrium_benchmark_se(x, par), x0, options);

theta_star = x(1);
sI_star    = x(2);
sU_star    = x(3);
w_star     = x(4);

q_star = qtheta(theta_star, mu, A);
f_star = ftheta(theta_star, mu, A);

% Hazards
fI_star = sI_star * f_star;
fU_star = sU_star * f_star;

% --- Steady-state shares (now hazards differ by state) ---
% i = delta*e/(fI + lambda)
% u = lambda*i/(fU)
% e+i+u=1
e_star = 1 / (1 + delta/(fI_star + lambda) + (lambda*delta)/(fU_star*(fI_star + lambda)));
i_star = delta * e_star / (fI_star + lambda);
u_star = lambda * i_star / fU_star;
u_total = i_star + u_star;

fprintf('\n---- Equilibrium Results (Endogenous Effort) ----\n');
fprintf('theta:        %.6f\n', theta_star);
fprintf('sI:           %.6f\n', sI_star);
fprintf('sU:           %.6f\n', sU_star);
fprintf('w:            %.6f\n', w_star);
fprintf('q(theta):     %.6f\n', q_star);
fprintf('f(theta):     %.6f\n', f_star);
fprintf('fI=sI*f:      %.6f\n', fI_star);
fprintf('fU=sU*f:      %.6f\n', fU_star);
fprintf('e:            %.6f\n', e_star);
fprintf('i:            %.6f\n', i_star);
fprintf('u:            %.6f\n', u_star);
fprintf('u_total:      %.6f\n', u_total);
fprintf('exitflag:     %d\n', exitflag);
fprintf('-----------------------------------------------\n');

%% SECTION 3: Comparative statics on bI (example)
par_old = par;

n = 10;
buffer_bI     = linspace(bI, bI + 0.10*bI, n);
buffer_theta  = NaN(1,n);
buffer_sI     = NaN(1,n);
buffer_sU     = NaN(1,n);
buffer_w      = NaN(1,n);
buffer_utotal = NaN(1,n);

options_cs = optimset('display','off', ...
                      'Algorithm','trust-region-dogleg', ...
                      'MaxIter',2000, ...
                      'Tolfun',1e-10, ...
                      'TolX',1e-10);

x0 = x; % warm start from equilibrium

for ii = 1:n
    par.bI = buffer_bI(ii);

    x = fsolve(@(x) equilibrium_benchmark_se(x, par), x0, options_cs);

    theta = x(1); sI = x(2); sU = x(3); w = x(4);

    q = qtheta(theta, mu, A);
    f = ftheta(theta, mu, A);

    fI = sI*f;
    fU = sU*f;

    e = 1 / (1 + delta/(fI + lambda) + (lambda*delta)/(fU*(fI + lambda)));
    i = delta * e / (fI + lambda);
    u = lambda * i / fU;

    buffer_theta(ii)  = theta;
    buffer_sI(ii)     = sI;
    buffer_sU(ii)     = sU;
    buffer_w(ii)      = w;
    buffer_utotal(ii) = i+u;

    x0 = x;
end

par = par_old;

%% Plots
figure(kpic); kpic=kpic+1;
plot(buffer_bI, buffer_theta, 'LineWidth', 1.5)
xlabel('Insured benefits (b_I)'); ylabel('Tightness (\theta)'); grid on; axis padded;

figure(kpic); kpic=kpic+1;
plot(buffer_bI, buffer_sI, 'LineWidth', 1.5)
xlabel('Insured benefits (b_I)'); ylabel('Search effort s_I'); grid on; axis padded;

figure(kpic); kpic=kpic+1;
plot(buffer_bI, buffer_sU, 'LineWidth', 1.5)
xlabel('Insured benefits (b_I)'); ylabel('Search effort s_U'); grid on; axis padded;

figure(kpic); kpic=kpic+1;
plot(buffer_bI, buffer_w, 'LineWidth', 1.5)
xlabel('Insured benefits (b_I)'); ylabel('Wage (w)'); grid on; axis padded;

figure(kpic); kpic=kpic+1;
plot(buffer_bI, buffer_utotal, 'LineWidth', 1.5)
xlabel('Insured benefits (b_I)'); ylabel('Total unemployment (i+u)'); grid on; axis padded;

% EOF
