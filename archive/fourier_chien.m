clear; close all; clc;

%% Fourier Curve
syms t
age = 20;
% Find parameterized x, y of Fourier Curve
[x_fourier, y_fourier] = fourier_curve(age);

% Obtain x,y coords
[X_fourier, Y_fourier] = fplot(x_fourier, y_fourier, [0, 2*pi]);
plot(X_fourier, Y_fourier)

%% Fitting Chien curve - Anterior
[X_fourierAnt, Y_fourierAnt] = fplot(x_fourier, y_fourier, [0, pi/2]);
figure
plot(X_fourierAnt, Y_fourierAnt)

% anterior a/b
a_ant = eval(subs(x_fourier, t, 0));
b0_ant = eval(subs(y_fourier, t, pi/2));

[b1_ant, b3_ant] = findChienCoefficients(X_fourierAnt(19:end), Y_fourierAnt(19:end), a_ant, b0_ant);

x_chienAnt = a_ant*sin(t);
y_chienAnt = (b0_ant + b1_ant*t^2 + b3_ant*t^4)*cos(t);

[X_chienAnt, Y_chienAnt] = fplot(x_chienAnt, y_chienAnt, [0 pi/2]);

hold on
plot(X_chienAnt, Y_chienAnt)


%% Ellipse - Anterior
x_elipAnt = a_ant*cos(t);
y_elipAnt = b0_ant*sin(t);

[X_elipAnt, Y_elipAnt] = fplot(x_elipAnt, y_elipAnt, [0 pi/2]);
plot(X_elipAnt, Y_elipAnt)

legend('Fourier', 'Chien', 'Ellipse')
title('Anterior Curve')

%% Anterior Energy (Fourier, Chien, Ellipse)
E_fourierAnt = findEnergy(x_fourier, y_fourier, pi/2, pi)

E_chienAnt = findEnergy(x_chienAnt, y_chienAnt, 0, pi/2)

E_elipAnt = findEnergy(x_elipAnt, y_elipAnt, 0, pi/2)

