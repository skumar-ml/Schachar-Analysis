clear; close all; clc;

%% Plot Urs curve
% For inner edge curve
urs_data_inner = readmatrix('Urs1_inner.csv');

% Take anterior portion & smooth
ursA_inner = urs_data_inner(urs_data_inner(:,1) < 0, :);
ursA_inner = sortrows(ursA_inner, 2);

% Split into X & Y (swapped so that we have a function)
Ya_inner = ursA_inner(:, 1)';
Xa_inner = ursA_inner(:, 2)';

% Smooth data & plot
Ya_inner_smoothed = smooth(Ya_inner)';
plot(Xa_inner, -Ya_inner_smoothed); hold on;
title("Anterior Curve (Urs 20 yr. old)")

% Repeat with outer edge curve
urs_data_outer = readmatrix('Urs1_outer.csv');

ursA_outer = urs_data_outer(urs_data_outer(:,1) < 0, :);
ursA_outer = sortrows(ursA_outer, 2);

Ya_outer = ursA_outer(:, 1)';
Xa_outer = ursA_outer(:, 2)';

Ya_outer_smoothed = smooth(Ya_outer)';
plot(Xa_outer, -Ya_outer_smoothed)

% Get final Urs curve by averaging outer & inner
X_Ant = (Xa_outer + Xa_inner)./2;
Y_Ant = (Ya_outer_smoothed + Ya_inner_smoothed)./2;

plot(X_Ant, -Y_Ant)
title("Raw data"); legend("Inner", "Outer", "Middle");


%% Forbes polynomial

% Center on origin
Y_Ant = Y_Ant - min(Y_Ant);

figure; plot(X_Ant, Y_Ant); hold on;

% Curvature of best-fit circle
syms t;

c = (2*max(Y_Ant) / (max(X_Ant)^2 + max(Y_Ant)^2));

r = 1/c;
x_circle = r*cos(t);
y_circle = r*sin(t)+r;

fplot(x_circle, y_circle, [-pi, 0]);
%% OLS Formulation
% Y
Y = Y_Ant'; % (d x 1)

% S
rho = X_Ant;
u = rho ./ max(rho);
M = 8;

lambda = (c * rho.^2) ./ (1+sqrt(1-c^2 .* rho.^2));
beta = (u.^2 .* (1 - u.^2)) ./ (sqrt(1-c^2.*rho.^2));

mini_S = (u.^2 .* ones(M+1, size(u,2))) .^ repmat((0:M)', 1, size(u,2));
mini_S = mini_S .* beta;

S = [lambda ; mini_S];

% A
A = inv(S * S') * S * Y;

%% Assembling curve 

% Points
Y_forbes = A'*S;
plot(X_Ant, Y_forbes)

% Equation (rho is equivalent to x, z is equivalent to y)
syms rho
u = rho / max(X_Ant);

z = (c*rho^2) / (1+sqrt(1-c^2*rho^2)) + (u^2 * (1-u^2)) / (sqrt(1-c^2*rho^2)) .* (A(2:end)' * (u.^2 .^ (0:M)'));
fplot(z, [-max(X_Ant), max(X_Ant)])

title("Best fit"); legend("Curve", "Best Fit Circle", "Forbes Polynomial (points)", "Forbes Polynomial (equation)");

% Difference
diff = abs(Y_forbes' - Y);

figure
plot(X_Ant, diff)
title("Sag difference")