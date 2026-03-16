%% Three-State Benchmark Search and Matching (SaM) model
%  Employed (E), Unemployed insured (I), Unemployed uninsured (U)
%  Transitions:
%    E -> I at rate delta
%    I -> U at rate lambda
%    I,U -> E at rate f(theta) = theta*q(theta)
%
%  Solves for (theta, w) using equilibrium_benchmark.m (three-state version)
%
% MScT DEPP (Franck Malherbet) — adapted to 3-state UI model

%% SECTION 1: Environment / Initialization / Parameters
clearvars; close all; clc;

% run from directory
cd(fileparts(mfilename('fullpath')))
addpath('Toolbox')

% Figure counter
kpic = 1;

%-----% Parameters / env. %-----------------------------------------------%
r_annual = 0.05;
r        = r_annual;

beta  = 0.5;    % bargaining weight
delta = 0.15;   % job destruction: E -> I
gamma = 0.15;   % vacancy posting cost parameter (as in benchmark)

% Productivity
y = 1;

% Matching function 
eta = 0.5;
m0  = 0.75;

% Policy parameters (three-state)
bI     = 0.25;  % insured unemployment benefit
bU     = 0.10;  % uninsured unemployment benefit (typically lower)
lambda = 0.50;  % expiration rate I -> U

% Create structure for parameters
par.name   = 'Parameters';
par.r      = r;
par.beta   = beta;
par.delta  = delta;
par.gamma  = gamma;
par.y      = y;
par.eta    = eta;
par.m0     = m0;

par.bI     = bI;
par.bU     = bU;
par.lambda = lambda;

%% SECTION 2: Solve for (aggregate) labor market tightness and the wage
options = optimset('display','iter', ...
                   'Algorithm','trust-region-dogleg', ...
                   'Diagnostics','off', ...
                   'MaxIter',1000, ...
                   'Tolfun',1e-8, ...
                   'TolX',1e-8);

x0 = [0.5 0.5];  % Initial guess: [theta, w]

[x, fval, exitflag] = fsolve(@(x) equilibrium_benchmark(x, par), x0, options);

theta_star = x(1);
w_star     = x(2);

% Matching objects at equilibrium
q_star  = qtheta(theta_star, eta, m0);     % vacancy-filling
f_star  = theta_star * q_star;             % job-finding

% --- Steady-state shares in 3-state model ---
% Flow balance:
% i = delta*e/(f+lambda)
% u = lambda*i/f
% e + i + u = 1
e_star = 1 / (1 + delta/(f_star + lambda) + (lambda*delta)/(f_star*(f_star + lambda)));
i_star = delta * e_star / (f_star + lambda);
u_star = lambda * i_star / f_star;

u_total = i_star + u_star;

% Display equilibrium results
fprintf('\n---- Equilibrium Results (3-state) ----\n');
fprintf('Market Tightness (theta):          %.6f\n', theta_star);
fprintf('Wage (w):                          %.6f\n', w_star);
fprintf('Vacancy-Filling Rate q(theta):     %.6f\n', q_star);
fprintf('Job-Finding Rate f(theta):         %.6f\n', f_star);
fprintf('Employment share e:                %.6f\n', e_star);
fprintf('Insured unemployment share i:      %.6f\n', i_star);
fprintf('Uninsured unemployment share u:    %.6f\n', u_star);
fprintf('Total unemployment i+u:            %.6f\n', u_total);
fprintf('Exitflag (fsolve):                 %d\n', exitflag);
fprintf('----------------------------------------\n');

%% SECTION 3: Comparative statics on insured benefits bI
% possible to play with U or lambda instead

% Save environment
par_old = par;

n = 10;
buffer_bI     = linspace(bI, bI + 0.10*bI, n);
buffer_theta  = NaN(1, n);
buffer_w      = NaN(1, n);
buffer_utotal = NaN(1, n);

options_cs = optimset('display','off', ...
                      'Algorithm','trust-region-dogleg', ...
                      'Diagnostics','off', ...
                      'MaxIter',1000, ...
                      'Tolfun',1e-8, ...
                      'TolX',1e-8);

x0 = [theta_star, w_star];  % warm start from equilibrium

for ii = 1:n
    par.bI = buffer_bI(ii);

    [x, ~, ~] = fsolve(@(x) equilibrium_benchmark(x, par), x0, options_cs);

    theta_cs = x(1);
    w_cs     = x(2);

    q_cs = qtheta(theta_cs, eta, m0);
    f_cs = theta_cs * q_cs;

    e_cs = 1 / (1 + delta/(f_cs + lambda) + (lambda*delta)/(f_cs*(f_cs + lambda)));
    i_cs = delta * e_cs / (f_cs + lambda);
    u_cs = lambda * i_cs / f_cs;

    buffer_theta(ii)  = theta_cs;
    buffer_w(ii)      = w_cs;
    buffer_utotal(ii) = i_cs + u_cs;

    x0 = x; % update warm start
end

% Restore environment
par = par_old;

%% Plots
figure(kpic); kpic = kpic + 1;
plot(buffer_bI, buffer_theta, 'LineWidth', 1.5)
xlabel('Insured benefits (b_I)')
ylabel('Labor market tightness (\theta)')
grid on; axis padded;

figure(kpic); kpic = kpic + 1;
plot(buffer_bI, buffer_w, 'LineWidth', 1.5)
xlabel('Insured benefits (b_I)')
ylabel('Wage (w)')
grid on; axis padded;

figure(kpic); kpic = kpic + 1;
plot(buffer_bI, buffer_utotal, 'LineWidth', 1.5)
xlabel('Insured benefits (b_I)')
ylabel('Total unemployment (i+u)')
grid on; axis padded;

% EOF
