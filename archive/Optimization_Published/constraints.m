function [c, ceq] = constraints( x )
%CONSTRAINTS Summary of this function goes here
%   Detailed explanation goes here

global divs;
global base;

% Check the arclength difference
c(1) = arc_diff(x) - .005;

% Ensure that the x axis length is 1.02 times the distance
ceq(1) = (1.02)*base(3) - x(3);

% Ensure symmetry at the axes
ceq(2) = x(1) - x(5);
ceq(3) = x(2) - x(4);
ceq(4) = x(divs) - x(divs-4);
ceq(5) = x(divs-1) - x(divs-3);

% Compute the necessary information for area constraint
area_b = sum_sectors(base);
area_x = sum_sectors(x);

% Ensure that the area is the same
ceq(6) = area_b - area_x;

end

