function [result] = kappa(r,r1,r2)
%KAPPA calculates the curvature based on an input radius
%      their derivatives at a specific angle

% Calculate the kappa function
numer = r^2 + 2*(r1^2) - (r*r2);
denom = (r^2 + r1^2)^(3/2);
result = numer/denom;

end

