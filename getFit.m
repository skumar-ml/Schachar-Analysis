function [fit] = getFit(raw_x, raw_y, eq)
% Finds fit of equation to raw data
syms t
syms y(t)

y(t) = eq;

eq_y = double(y(raw_x)); % Value of equation evaluated at all points

% Calculate root mean-sqaured fit

fit = sqrt((mean((eq_y - raw_y)).^2) * 10^6);
end

