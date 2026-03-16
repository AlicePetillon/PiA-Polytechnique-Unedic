% file: Toolbox/ftheta.m
function f = ftheta(theta, mu, A)
% f(theta) = A * theta^(1-mu)
f = A .* (theta .^ (1 - mu));
end
