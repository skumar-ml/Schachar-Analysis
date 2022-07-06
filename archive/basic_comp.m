clear; clc; close all;
%% Excel Curve
syms t
x_excel = t;
y_excel = -0.0022*t^6 + 0.0251*t^5 - 0.1074*t^4 + 0.21*t^3 - 0.2446*t^2 + 0.065*t + 1.5926;

x_p = diff(x_excel, t, 1);
x_pp = diff(x_excel, t, 2);
y_p = diff(y_excel, t, 1);
y_pp = diff(y_excel, t, 2);

k_excel = (x_p*y_pp - y_p*x_pp) / (x_p^2 + y_p^2)^(3/2);

expr_excel = k_excel^2 * ((x_p^2 + y_p^2) ^ (1/2));
E_excel = vpaintegral(expr_excel, 0, 4.368)

figure
[X_excel, Y_excel] = fplot(x_excel, y_excel, [0 4.368]);
plot(X_excel, Y_excel)

% PART 2 %
syms t
x_excel = 4.368*sin(t);
y_excel = -0.0022*(x_excel)^6 + 0.0251*(x_excel)^5 - 0.1074*(x_excel)^4 + 0.21*(x_excel)^3 - 0.2446*(x_excel)^2 + 0.065*(x_excel) + 1.5926;

x_p = diff(x_excel, t, 1);
x_pp = diff(x_excel, t, 2);
y_p = diff(y_excel, t, 1);
y_pp = diff(y_excel, t, 2);

k_excel = (x_p*y_pp - y_p*x_pp) / (x_p^2 + y_p^2)^(3/2);

expr_excel = k_excel^2 * ((x_p^2 + y_p^2) ^ (1/2));
E_excel = vpaintegral(expr_excel, 0, pi/2)

hold on
[X_excel, Y_excel] = fplot(x_excel, y_excel, [0 pi/2]);
plot(X_excel, Y_excel)

%% Dr. Ron's Curve
syms t
b0 = 1.5953805;
b1 = -0.0962011;
b3 = 0.11256773;
a = 4.3682554;

x_ron =  a*sin(t);
y_ron = (b0 + b1*t^2 + b3*t^4)*cos(t);

% pre-calculating first/second order derivatives of x(t), y(t)
tic
x_p = diff(x_ron, t, 1);
x_pp = diff(x_ron, t, 2);
y_p = diff(y_ron, t, 1);
y_pp = diff(y_ron, t, 2);

% formula for a curvature of a curve parameterized by x(t), y(t) -
% wikipedia

k_ron = (x_p*y_pp - y_p*x_pp) / (x_p^2 + y_p^2)^(3/2);

% expression to be integrated -- from Horn pg. 5
expr_ron = k_ron^2 * ((x_p^2 + y_p^2) ^ (1/2));

E_ron = vpaintegral(expr_ron, 0, pi/2)

hold on
[X_ron, Y_ron] = fplot(x_ron, y_ron, [0 pi/2]);
plot(X_ron, Y_ron)
% 
% [K_ron, E_ron] = findEnergy(X_ron, Y_ron);
% figure
% ezplot(K_ron, [min(X_ron) max(X_ron)])
% title("Dr. Ron Curvature")

%% Optimal Ellipse 
% major/minor axis of the ellipse
% Gen rule - a = a, b = b0
% to verify Horn's results, set a=1, b=1.32
% for Strenk Anteroir A, set a=4.5, b=1.44
% a = 4.5; 
b = b0;

% parameterization of x(t), y(t) -- from Horn pg. 4
x_elip = a*cos(t);
y_elip = b*sin(t);

eccentricity = sqrt(1 - (a/b)^2);

k_elip = (a*b) / (a^2*(sin(t))^2 + b^2*(cos(t))^2)^(3/2);

% expression to be integrated -- Horn pg. 5
expr = (1 / (1 - eccentricity^2*(sin(t))^2)^(5/2));
% expr_2 = k_elip^2 * sqrt(diff(x_elip,t,1)^2 + diff(y_elip,t,1)^2);

E_elip = (a^2/b^3) * vpaintegral(expr, 0, pi/2)
% E_elip2 = vpa(int(expr_2, 0, pi/2))
toc

[X_elip, Y_elip] = fplot(x_elip, y_elip, [0, pi/2]);
plot(X_elip, Y_elip)
title("Curves | Lens 22year")
legend("Excel", "Excel pt. 2", "Chien", "Ellipse")
