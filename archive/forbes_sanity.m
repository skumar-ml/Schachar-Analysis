close all; clear all; clc;

%% Define test data
syms x;
y = 1/40 * x^2;
figure(1); fplot(y, [-20, 20])
[X_data, Y_data] = fplot(y, [-20, 20], 'MeshDensity', 100); hold on;

%% My Urs Implementation

% Find forbes equation
[forbes_eq, Y_forbes, A_mine] = forbes(X_data, Y_data, 6);
A_mine = A_mine';
plot(X_data, Y_forbes); 

% Fit Error
diff = Y_forbes - Y_data;
figure(2); plot(X_data, diff)
xlim([0, 20])

% Mean Squared Difference
fit_myUrs = mean((Y_forbes - Y_data)).^2; 

%% Paper Urs Implementation
% Coefficients
A_urs = [2019004, 7143, -13944, 4190, -1095, 283, -68] * 10^-9; % from Robust paper

% Set up equation
syms rho;
u = rho / max(X_data);
c = (2*max(Y_data) / (max(X_data)^2 + max(Y_data)^2)); 
M = size(A_urs,2)-1;

z = (c*rho^2) / (1+sqrt(1-c^2*rho^2)) + (u^2 * (1-u^2)) / (sqrt(1-c^2*rho^2)) .* (A_urs * (u.^2 .^ ((0:M)')));

% Evaluate at x data points
FUN = matlabFunction(z);
Y_forbesPaper = feval(FUN, X_data);
figure(1); plot(X_data, Y_forbesPaper); legend("Raw", "Fourier (Mine)", "Fourier (Paper)")

% Fit Error
diff = abs(Y_forbesPaper - Y_data);
figure(3); plot(X_data, diff)
xlim([0, 20])

% Mean S
% Mean Squared Difference
fit_paperUrs = mean((Y_forbesPaper - Y_data)).^2; 




