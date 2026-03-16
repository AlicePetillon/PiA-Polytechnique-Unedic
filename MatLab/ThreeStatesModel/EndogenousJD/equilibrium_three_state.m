function System = equilibrium_three_state(x, par)


theta = x(1);
R     = x(2);

% Parameters 
r       = par.r;
mu      = par.mu;       % Shock frequency 
lambda  = par.lambda;   % Transition I -> U 
beta    = par.beta;
gamma   = par.gamma;    % Vacancy cost 
yub     = par.yub;      
bI      = par.bI;
bU      = par.bU;

% Intermediate Matching Functions
q_theta = par.m0 * theta^(-par.eta); 

% Surplus Integral (Assuming Uniform Distribution)
surplus_integral = integral(@(y) (y - R) .* unifpdf(y, par.ylb, yub), R, yub);

% Value of Insured Unemployment 
term_VU_VI = (lambda * (bU - bI)) / (r + lambda + (theta * gamma * (r + mu)) / ((1 - beta) * (yub - R)));
rVI = bI + (theta * beta * gamma) / (1 - beta) + term_VU_VI;

% System of Equations
System = [   
    % EQ 1: Job Creation 
    (gamma / q_theta) - ((1 - beta) * (yub - R) / (r + mu));

    % EQ 2: Job Destruction 
    R - rVI + (mu / (r + mu)) * (surplus_integral / (1)); 
];
end