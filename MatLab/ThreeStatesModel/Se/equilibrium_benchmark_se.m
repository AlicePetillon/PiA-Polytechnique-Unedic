% file: equilibrium_benchmark_se.m
function System = equilibrium_benchmark_se(x, par)
% Three-State DMP with Time-Limited UI and Endogenous Search Effort
% Unknowns: theta, sI, sU, w
% Equilibrium system

theta = x(1);
sI    = x(2);
sU    = x(3);
w     = x(4);

% Parameters
r      = par.r;
beta   = par.beta;
delta  = par.delta;
kappa  = par.kappa;
y      = par.y;

mu     = par.mu;
A      = par.A;

bI     = par.bI;
bU     = par.bU;
lambda = par.lambda;

c0     = par.c0;

% Matching
q = qtheta(theta, mu, A);
f = ftheta(theta, mu, A);

% Nash surplus piece = (beta/(1-beta)) * kappa/q
S = (beta/(1 - beta)) * (kappa / q);

% (i) Job creation: kappa/q = (y-w)/(r+delta)
eq_JC = kappa / q - (y - w) / (r + delta);

% (ii) Insured effort: sI = (f/c0)*S
eq_sI = sI - (f / c0) * S;

% (iii) Uninsured effort quadratic (eq 20)
eq_sU = 0.5*c0*sU^2 + c0*(r + lambda)/f * sU + ...
        ( (bU - bI) - (r + lambda)*S - 0.5*(f^2/c0)*S^2 );

% (iv) Wage equation (eq 21)
w_rhs = bI - 0.5*c0*sI^2 + (r + delta + sI*f + lambda)*S - lambda*(c0/f)*sU;
eq_W  = w - w_rhs;

System = [eq_JC; eq_sI; eq_sU; eq_W];

end
