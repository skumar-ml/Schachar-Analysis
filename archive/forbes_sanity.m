close all; clear all; clc;

% Define test model
syms x;
y = 1/40 * x^2;
fplot(y, [-20, 20])
[X_data, Y_data] = fplot(y, [-20, 20], 'MeshDensity', 100); hold on;

% Find forbes equation
[forbes_eq, Y_forbes, A] = forbes(X_data, Y_data, 8)
plot(X_data, Y_forbes); legend("Raw", "Fourier")

% Difference
diff = abs(Y_forbes - Y_data);
figure; plot(X_data, diff)
xlim([0, 20])


