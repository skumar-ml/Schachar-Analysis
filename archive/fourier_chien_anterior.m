clear; close all;

%% Fourier Curve
syms t; % t is in radians.
age = 20;
% Find parameterized x, y of Fourier Curve
[x_fourier, y_fourier] = fourier_curve(age);

% Obtain x,y coords
[Y_fourier, X_fourier] = fplot(x_fourier, y_fourier, [0, 2*pi]); % swapped to make anterior and posterior an actual function
plot(X_fourier, -Y_fourier) % flat = anterior, curved = posterior | negative flips orientation of lens
title("Entire Fourier Lens (Flat = anterior)")
%% Fitting Chien curve - Anterior
fourier_bounds = [pi/2, 3*pi/2];
[Y_fourierAnt, X_fourierAnt] = fplot(x_fourier, y_fourier, fourier_bounds);
figure
plot(X_fourierAnt, -Y_fourierAnt)

% anterior a/b
b0_ant = min(Y_fourierAnt);
a_ant = max(X_fourierAnt);

[b1_ant, b3_ant] = findChienCoefficients(X_fourierAnt, Y_fourierAnt, a_ant, b0_ant);

x_chienAnt = a_ant*sin(t); % in mm
y_chienAnt = (b0_ant + b1_ant*t^2 + b3_ant*t^4)*cos(t); % in mm

chien_bounds = [-pi/2, pi/2];
[X_chienAnt, Y_chienAnt] = fplot(x_chienAnt, y_chienAnt, chien_bounds);

hold on
plot(X_chienAnt, -Y_chienAnt)


%% Ellipse - Anterior
x_elipAnt = a_ant*cos(t); % in mm
y_elipAnt = b0_ant*sin(t); % in mm

elip_bounds = [0 pi];
[X_elipAnt, Y_elipAnt] = fplot(x_elipAnt, y_elipAnt, elip_bounds);
plot(X_elipAnt, -Y_elipAnt)

legend('Fourier', 'Chien', 'Ellipse')
title('Anterior Curve')
ylabel('Lens Thickness (mm)')
xlabel('Distance from Optic Axis (mm)')

%% Anterior Energy & Variation (Fourier, Chien, Ellipse)
% offset = 0.5;
zone = 3; % optical zone (in mm)
offset = fourier_bounds(1) - asin(zone/a_ant); % Convert from mm to radians

[E_fourierAnt, k_fourierAnt] = findEnergy(x_fourier, y_fourier, fourier_bounds(1)+offset, fourier_bounds(2)-offset)
[E_chienAnt, k_chienAnt] = findEnergy(x_chienAnt, y_chienAnt, chien_bounds(1)+offset, chien_bounds(2)-offset)
[E_elipAnt, k_elipAnt] = findEnergy(x_elipAnt, y_elipAnt, elip_bounds(1)+offset, elip_bounds(2)-offset)

[bendE_elipAnt, bendExpr_elipAnt] = findBendingEnergy(x_elipAnt, y_elipAnt, elip_bounds(1)+offset, elip_bounds(2)-offset);
figure
fplot(x_elipAnt, bendExpr_elipAnt, [elip_bounds(1)+offset, elip_bounds(2)-offset])
% findBendingEnergy(x_chienAnt, y_chienAnt, chien_bounds(1)+offset, chien_bounds(2)-offset)
% findBendingEnergy(x_fourier, y_fourier, fourier_bounds(1)+offset, fourier_bounds(2)-offset)

% Graph of curvatures - x-axis in t (converted from radians)
figure
fplot(a_ant*sin(t - pi/2 - pi/2), abs(k_fourierAnt), fourier_bounds)
hold on
fplot(a_ant*sin(t + pi/2 - pi/2), abs(k_chienAnt), chien_bounds)
fplot(a_ant*sin(t - pi/2), abs(k_elipAnt), elip_bounds)

legend("Fourier", "Chien", "Ellipse")
title("Radius of Curvature for anterior")
% ylabel("Lens Curvature (mm^-1)")
ylabel("Radius of Curvature (mm)")
ylim([0, 20])
xlabel('Distance from Optic Axis (mm)')

