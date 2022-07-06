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
title("Anterior Curve (Urs 20 yr. old)")

% Get final Urs curve by averaging outer & inner
X_Ant = (Xa_outer + Xa_inner)./2;
Y_Ant = (Ya_outer_smoothed + Ya_inner_smoothed)./2;
plot(X_Ant, -Y_Ant)

%% Plot Fourier Curve
age = 20;
% Find parameterized x, y of Fourier Curve
[x_fourier, y_fourier] = fourier_curve(age);

% Obtain x,y coords
[X_fourier, Y_fourier] = fplot(x_fourier, y_fourier, [0, 2*pi]); % swapped to make anterior and posterior an actual function
plot(X_fourier, -Y_fourier)
figure

%% Fitting Chien curve - Anterior
syms t
fourier_bounds = [pi/2, 3*pi/2];
[X_fourierAnt, Y_fourierAnt] = fplot(x_fourier, y_fourier, fourier_bounds);
plot(X_fourierAnt, -Y_fourierAnt)
legend("Raw Outer", "Raw Inner", "Raw Middle", "Fourier Curve"); xlabel("mm"); ylabel("mm")

% Fit chien to raw
b0_ant = min(Y_Ant);
a_ant = max(X_Ant) + 0.0001; % add epsilon for numerical stability

[b1_ant, b3_ant] = findChienCoefficients(X_Ant, Y_Ant, a_ant, b0_ant);

x_chienAnt = a_ant*sin(t);
y_chienAnt = (b0_ant + b1_ant*t^2 + b3_ant*t^4)*cos(t);

chien_bounds = [-pi/2, pi/2];
[X_chienAnt, Y_chienAnt] = fplot(x_chienAnt, y_chienAnt, chien_bounds);

figure
hold on
plot(X_chienAnt, -Y_chienAnt)
plot(X_fourierAnt, -Y_fourierAnt)
plot(X_Ant, -Y_Ant)

% Fit chien to fourier
b0_ant = min(Y_fourierAnt);
a_ant = max(X_fourierAnt) + 0.0001; % add epsilon for numerical stability

[b1_ant, b3_ant] = findChienCoefficients(X_fourierAnt, Y_fourierAnt, a_ant, b0_ant);

x_chienAnt2 = a_ant*sin(t);
y_chienAnt2 = (b0_ant + b1_ant*t^2 + b3_ant*t^4)*cos(t);

chien_bounds = [-pi/2, pi/2];
[X_chienAnt2, Y_chienAnt2] = fplot(x_chienAnt2, y_chienAnt2, chien_bounds);
plot(X_chienAnt2, -Y_chienAnt2)

%% Ellipse - Anterior
x_elipAnt = a_ant*cos(t); % in mm
y_elipAnt = b0_ant*sin(t); % in mm

elip_bounds = [0 pi];
[X_elipAnt, Y_elipAnt] = fplot(x_elipAnt, y_elipAnt, elip_bounds);
plot(X_elipAnt, -Y_elipAnt)

legend("Chien (fit to raw)", "Fourier", "Raw Data", "Chien (fit to fourier)", "Ellipse"); xlabel("mm"); ylabel("mm")

%% Anterior Energy & Variation (Fourier, Chien, Ellipse)
% offset = 0.5;
zone = 3; % optical zone [-zone, +zone]
offset = -1*(asin(zone/a_ant) - fourier_bounds(1)); % convert to polar

[E_fourierAnt, k_fourierAnt] = findEnergy(x_fourier, -y_fourier, fourier_bounds(1)+offset, fourier_bounds(2)-offset)
[E_elipAnt, k_elipAnt] = findEnergy(x_elipAnt, -y_elipAnt, elip_bounds(1)+offset, elip_bounds(2)-offset)
[E_chienAnt, k_chienAnt] = findEnergy(x_chienAnt, -y_chienAnt, chien_bounds(1)+offset, chien_bounds(2)-offset)
[E_chienAnt2, k_chienAnt2] = findEnergy(x_chienAnt2, -y_chienAnt2, chien_bounds(1)+offset, chien_bounds(2)-offset)

smth_fourierAnt = eval(vpaintegral(diff(k_fourierAnt, t, 2) ^ 2, fourier_bounds(1)+offset, fourier_bounds(2)-offset));
smth_elipAnt = eval(vpaintegral(diff(k_elipAnt, t, 2) ^ 2, elip_bounds(1)+offset, elip_bounds(2)-offset));
smth_chienAnt = eval(vpaintegral(diff(k_chienAnt, t, 2) ^ 2, chien_bounds(1)+offset, chien_bounds(2)-offset));
smth_chienAnt2 = eval(vpaintegral(diff(k_chienAnt2, t, 2) ^ 2, chien_bounds(1)+offset, chien_bounds(2)-offset));

figure; hold on;
fplot(abs(k_fourierAnt), fourier_bounds)
fplot(abs(diff(k_fourierAnt, t, 1)), fourier_bounds)

[bendE_chienAnt, firstD_chienAnt, expr_chienAnt] = findBendingEnergy(x_chienAnt, -y_chienAnt, chien_bounds(1)+offset, chien_bounds(2)-offset)
% figure; fplot(expr_chienAnt, [chien_bounds]); title("Chien Ant")
% 
[bendE_chienAnt2, firstD_chienAnt2, expr_chienAnt2] = findBendingEnergy(x_chienAnt2, -y_chienAnt2, chien_bounds(1)+offset, chien_bounds(2)-offset)
% figure; fplot(expr_chienAnt2, [chien_bounds]); title("Chien Ant 2")
% 
[bendE_elipAnt, firstD_elipAnt, expr_elipAnt] = findBendingEnergy(x_elipAnt, -y_elipAnt, elip_bounds(1)+offset, elip_bounds(2)-offset)
% figure; fplot(expr_elipAnt, [elip_bounds]); title("Ellipse Ant")
% 
[bendE_fourierAnt, firstD_fourierAnt, expr_fourierAnt] = findBendingEnergy(x_fourier, -y_fourier, fourier_bounds(1)+offset, fourier_bounds(2)-offset)
% figure; fplot(expr_fourierAnt, [fourier_bounds]); title("Fourier Ant")

% Graph of curvatures
figure; hold on;
fplot(a_ant*sin(t - pi/2 - pi/2), abs(k_fourierAnt), fourier_bounds)
fplot(a_ant*sin(t - pi/2), abs(k_elipAnt), elip_bounds)
fplot(a_ant*sin(t + pi/2 - pi/2), abs(k_chienAnt), chien_bounds)
fplot(a_ant*sin(t + pi/2 - pi/2), abs(k_chienAnt2), chien_bounds)

%avg_roc_fourier = eval(vpaintegral(1/k_fourierAnt, t, fourier_bounds(1)+offset, fourier_bounds(2)-offset)) / abs(fourier_bounds(1)+offset - fourier_bounds(2)-offset)
%avg_roc_chien = eval(vpaintegral(1/k_chienAnt, t, chien_bounds(1)+offset, chien_bounds(2)-offset)) / abs(chien_bounds(1)+offset - chien_bounds(2)-offset)

title("Radius of Curvature"); legend("Fourier", "Ellipse", "Chien (Raw)", "Chien (Fourier)"); xlabel("Distance from Optic Axis (mm)");
% title("ds (arc length)"); legend("Fourier", "Ellipse", "Chien (Raw)", "Chien (Fourier)"); xlabel("Distance from Optic Axis (mm)"); ylabel("Arc Length");

%% Plots

figure; hold on;
fplot(x_fourier, -y_fourier, fourier_bounds)
fplot(x_chienAnt, -y_chienAnt, chien_bounds)
legend("fourier", "chien")
