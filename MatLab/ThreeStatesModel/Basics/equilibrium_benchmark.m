function System = equilibrium_benchmark(x, par)
% Three-state DMP with time-limited UI
% Unknowns: theta, w
% Equations:
% (1) Wage setting
% (2) Job creation: gamma/q(theta) = (y - w)/(r + delta)

%% Endogenous variables
theta = x(1);
w     = x(2);

%% Parameters
r      = par.r;
beta   = par.beta;
delta  = par.delta;
gamma  = par.gamma;
y      = par.y;

eta    = par.eta;
m0     = par.m0;

% Three-state policy parameters
bI     = par.bI;       % insured unemployment benefit
bU     = par.bU;       % uninsured unemployment benefit
lambda = par.lambda;   % expiration rate I -> U

%% Matching objects
q = qtheta(theta, eta, m0);   % vacancy filling rate q(theta)
f = theta * q;                % job finding rate f(theta) = theta*q(theta)

%% Wage-setting curve: WSfinal
% B(theta) = (bU - bI)/(r + f(theta) + lambda)
B = (bU - bI) / (r + f + lambda);

% Numerator:
% beta*y + (1-beta)*(r+delta)/(r+delta+f) * [bI + lambda*B]
num = beta * y + (1 - beta) * (r + delta) / (r + delta + f) * (bI + lambda * B);

% Denominator:
% 1 - (1-beta)*f/(r+delta+f)
den = 1 - (1 - beta) * f / (r + delta + f);

w_theta = num / den;

%% System of equilibrium conditions
System = [
    % (1) Wage equation
    w - w_theta;

    % (2) Labor market tightness / job creation
    gamma / q - (y - w) / (r + delta)
];

end
