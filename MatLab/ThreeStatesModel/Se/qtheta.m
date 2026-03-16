% file: Toolbox/qtheta.m
function q = qtheta(theta, mu, A)
% q(theta) = A * theta^(-mu)
q = A .* (theta .^ (-mu));
end
