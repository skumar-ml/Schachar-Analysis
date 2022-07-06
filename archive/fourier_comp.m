clear; clc; close all;
%% Fourier Curve
% Defining coefficients
syms t
A_n1 = [2.6466 0.2246 -0.97938 0.010573 0.37993 -0.032321 -0.16846 0.027934 0.066522 -0.014232 -0.021375];
A_n2 = [812.11 170.62 -297.37 -34.901 -26.276 1.6647 69.192 -9.5571 -42.251 1.7295 18.638] .* 1e-5;
age = 20;

A = A_n1 + A_n2 .* age;

% Defining equation
rho = t-t;

for i=0:10
    rho = rho + A(i+1)*cos(i*t);
end

% curve is polar, paramaterized in terms of cartesian

% Parameterize the curve & plot
x_fourier = rho * cos(t);
y_fourier = rho * sin(t);

[X_fourier, Y_fourier] = fplot(x_fourier, y_fourier, [0, 2*pi]);

% Pre-calculate derivatives
x_p = diff(x_fourier, t, 1);
x_pp = diff(x_fourier, t, 2);
y_p = diff(y_fourier, t, 1);
y_pp = diff(y_fourier, t, 2);

k_fourier = (x_p*y_pp - y_p*x_pp) / (x_p^2 + y_p^2)^(3/2);

expr_fourier = k_fourier^2 * ((x_p^2 + y_p^2) ^ (1/2));

E_fourier = vpaintegral(expr_fourier, pi/2, 3*pi/2)

%% Ellipse Curve
syms t
a = subs(x_fourier, t, pi); % x-axis
b = subs(y_fourier, t, pi/2); % y-axis

% parameterization of x(t), y(t) -- from Horn pg. 4
x_elip = a*cos(t);
y_elip = b*sin(t);

% Pre-calculate derivatives
x_p = diff(x_elip, t, 1);
x_pp = diff(x_elip, t, 2);
y_p = diff(y_elip, t, 1);
y_pp = diff(y_elip, t, 2);

k_elip = (x_p*y_pp - y_p*x_pp) / (x_p^2 + y_p^2)^(3/2);


expr_elip = k_elip^2;

E_elip = vpaintegral(expr_elip, -pi/2, pi/2)

[X_elip, Y_elip] = fplot(x_elip, y_elip, [0, 2*pi]);


%% Plot of fourier & ellipse

% Plot of curve

figure
plot(X_fourier, Y_fourier)
hold on
plot(X_elip, Y_elip)
legend("Fouier", "Ellipse")
title("Lens")

% Plot of expression that is integrated
figure
fplot(t, expr_fourier, [0, 2*pi])
hold on
fplot(t, expr_elip, [0, 2*pi])
legend("Fouier", "Ellipse")
title("Derivative of energy (expression that is integrated)")

% Plot of curvature
figure
fplot(t, k_fourier.^2, [0, 2*pi])
hold on
fplot(t, k_elip.^2, [0, 2*pi])
legend("Fouier", "Ellipse")
title("Curvature squared (kappa)")







