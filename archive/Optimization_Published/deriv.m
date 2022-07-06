function [ d ] = deriv(r, i, j)
%DERIV Take the derivatives of r and generate
%      the results in an array

global divs;
global theta;

% Initialize the set for the derivative
d = zeros(1, divs);

% Calculate the derivative for each angle specified
for k=i:j

    r1 = r(k-1); 
    r2 = r(k+1);

    
    t1 = theta(k-1);
    t2 = theta(k+1);

    slope1 = (r2 - r1)/(t2 - t1);

    % Derivative is delta radius / delta angle
    d(k) = slope1;

end

end
