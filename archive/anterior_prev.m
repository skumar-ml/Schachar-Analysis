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

%% Obtain & Plot Fourier Curve
age = 20;
% Find parameterized x, y of Fourier Curve
[x_fourier, y_fourier] = fourier_curve(age);

% Obtain x,y coords
[X_fourier, Y_fourier] = fplot(x_fourier, y_fourier, [0, 2*pi]); % swapped to make anterior and posterior an actual function
plot(X_fourier, -Y_fourier)
figure

% Extract anterior portion of Fourier curve
fourier_bounds = [pi/2, 3*pi/2];
[X_fourierAnt, Y_fourierAnt] = fplot(x_fourier, y_fourier, fourier_bounds);
plot(X_fourierAnt, -Y_fourierAnt)
legend("Raw Outer", "Raw Inner", "Raw Middle", "Fourier Curve"); xlabel("mm"); ylabel("mm")

%% Fitting Chien curve - Anterior
syms t

% Fit chien to raw
b0_ant = min(Y_Ant);
a_ant = max(X_Ant) + 0.0001; % add epsilon for numerical stability

[b1_ant, b3_ant] = findChienCoefficients(X_Ant, Y_Ant, a_ant, b0_ant);

x_chienAnt = a_ant*sin(t);
y_chienAnt = (b0_ant + b1_ant*t^2 + b3_ant*t^4)*cos(t);

chien_bounds = [-pi/2, pi/2];
[X_chienAnt, Y_chienAnt] = fplot(x_chienAnt, y_chienAnt, chien_bounds);

%% Fitting Forbes curve - Anterior
syms rho;

Y_center = Y_Ant - min(Y_Ant);
figure; plot(X_Ant, Y_center); title("Centered raw data curve"); hold  on;

[forbes_eq_centered, Y_forbes_centered] = forbes(X_Ant, Y_center, 3);

Y_forbes = -1 * (Y_forbes_centered + min(Y_Ant));

forbes_eq = -forbes_eq_centered;
forbes_eq = subs(forbes_eq, rho, t); % done because functions use t as sym variable
% Note -- the t in the forbes equation is cartesian! Stands for x (not
% theta)
%% Ellipse - Anterior
x_elipAnt = a_ant*cos(t); % in mm
y_elipAnt = b0_ant*sin(t); % in mm

elip_bounds = [0 pi];
[X_elipAnt, Y_elipAnt] = fplot(x_elipAnt, y_elipAnt, elip_bounds);

% Plot raw data, chien, fourier, & ellipse anterior curves
figure; hold on
plot(X_chienAnt, -Y_chienAnt); plot(X_fourierAnt, -Y_fourierAnt); plot(X_Ant, -Y_Ant); plot(X_elipAnt, -Y_elipAnt); plot(X_Ant, Y_forbes);
legend("Chien (fit to raw)", "Fourier", "Raw Data", "Ellipse", "Forbes"); xlabel("mm"); ylabel("mm")

%% Anterior Energy & Variation (Fourier, Chien, Ellipse)
zone = 3; % optical zone = zone * 2 (in mm)
offset = -1*(asin(zone/a_ant) - fourier_bounds(1)); % convert to polar

% Find curvature
k_chienAnt = findCurvature(x_chienAnt, -y_chienAnt, chien_bounds(1)+offset, chien_bounds(2)-offset);
k_elipAnt = findCurvature(x_elipAnt, -y_elipAnt, elip_bounds(1)+offset, elip_bounds(2)-offset);
k_fourierAnt = findCurvature(x_fourier, -y_fourier, fourier_bounds(1)+offset, fourier_bounds(2)-offset);
k_forbes = findCurvature(t, forbes_eq, -zone, zone);

% Find smoothing energy (integral of derivative of curvature squared)
smth_chienAnt = eval(vpaintegral(diff(k_chienAnt, t, 1) ^ 2, chien_bounds(1)+offset, chien_bounds(2)-offset));
smth_elipAnt = eval(vpaintegral(diff(k_elipAnt, t, 1) ^ 2, elip_bounds(1)+offset, elip_bounds(2)-offset));
smth_fourierAnt = eval(vpaintegral(diff(k_fourierAnt, t, 1) ^ 2, fourier_bounds(1)+offset, fourier_bounds(2)-offset));
smth_forbes = eval(vpaintegral(diff(k_forbes, t, 1) ^ 2, -zone, zone));

% Find bending energy
%[bendE_chienAnt, firstD_chienAnt, expr_chienAnt] = findBendingEnergy(x_chienAnt, -y_chienAnt, chien_bounds(1)+offset, chien_bounds(2)-offset);
[bendE_elipAnt, firstD_elipAnt, expr_elipAnt] = findBendingEnergy(x_elipAnt, -y_elipAnt, elip_bounds(1)+offset, elip_bounds(2)-offset);
%[bendE_fourierAnt, firstD_fourierAnt, expr_fourierAnt] = findBendingEnergy(x_fourier, -y_fourier, fourier_bounds(1)+offset, fourier_bounds(2)-offset);
%[bendE_forbes, firstD_forbes, expr_forbes] = findBendingEnergy(t, forbes_eq, -zone, zone);

% Mean of RoC
meanROC_chienAnt = abs(1/((chien_bounds(2)-offset) - (chien_bounds(1)+offset)) * eval(vpaintegral(1/k_chienAnt, chien_bounds(1)+offset, chien_bounds(2)-offset)));
valROC_chienAnt = abs(eval(subs(1/k_chienAnt, t, chien_bounds(1)+offset)));

meanROC_elipAnt = abs(1/((elip_bounds(2)-offset) - (elip_bounds(1)+offset)) * eval(vpaintegral(1/k_elipAnt, elip_bounds(1)+offset, elip_bounds(2)-offset)));
valROC_elipAnt = abs(eval(subs(1/k_elipAnt, t, elip_bounds(1)+offset)));

meanROC_fourierAnt = abs(1/((fourier_bounds(2)-offset) - (fourier_bounds(1)+offset)) * eval(vpaintegral(1/k_fourierAnt, fourier_bounds(1)+offset, fourier_bounds(2)-offset)));
valROC_fourierAnt = abs(eval(subs(1/k_fourierAnt, t, fourier_bounds(1)+offset)));

%% Graphs
% Graph of derivative of curvature
figure; hold on; fplot(abs(diff(k_chienAnt, t, 1)), chien_bounds); fplot(abs(diff(k_fourierAnt, t, 1)), fourier_bounds);
title("Absolute Value of the Derivative of curvature"); legend("Chien", "Fourier");

% Graph of curvatures
figure; hold on;
fplot(a_ant*sin(t - pi/2 - pi/2), abs(1/k_fourierAnt), fourier_bounds)
fplot(a_ant*sin(t - pi/2), abs(1/k_elipAnt), elip_bounds)
fplot(a_ant*sin(t + pi/2 - pi/2), abs(1/k_chienAnt), chien_bounds)
fplot(t, abs(1/k_forbes))

title("Radius of Curvature"); legend("Fourier", "Ellipse", "Chien (Raw)", "Forbes"); xlabel("Distance from Optic Axis (mm)");