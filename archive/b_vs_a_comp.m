clear; clc; close all;
%% Comparison
ratio = 0.05 : 0.05 : 1; % ratio of b0/a

E_ron_total = zeros(1,length(ratio));
E_elip_total = zeros(1, length(ratio));

for i=1:length(ratio)
    i
    % Define parameters - rest of numbers taken from Strenk A, Anterior
    a = 4.5;
    
    b0 = a*ratio(i);
    b1 = -0.0962011;
    b3 = 0.11256773;
 
    % Calculate energy for Chien and Ellipse
    [E_ron, x_ron, y_ron] = ron_curve(a, b0, b1, b3);
    [E_elip, x_elip, y_elip] = elip_curve(a, b0);
    
    % Save
    E_ron_total(i) = E_ron; 
    E_elip_total(i) = E_elip;

    % plot
    figure
    [X_ron, Y_ron] = fplot(x_ron, y_ron, [0 pi/2]);
    plot(X_ron, Y_ron)
    hold on
    [X_elip, Y_elip] = fplot(x_elip, y_elip, [0, pi/2]);
    plot(X_elip, Y_elip)
    title(strcat("Strenk A, Anterior | b0: ", int2str(b0)))
    legend("Chien", "Ellipse")
end

save('strenkA_anterior.mat', 'E_ron_total', 'E_elip_total')

figure
hold on
plot(ratio, E_ron_total)
plot(ratio, E_elip_total)

%% Functions
function [E_ron, x_ron, y_ron] = ron_curve(a, b0, b1, b3)
syms t
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
expr_ron = k_ron^2 * (x_p^2 + y_p^2) ^ (1/2);

E_ron = vpa(int(expr_ron, 0, pi/2));
end
function [E_elip, x_elip, y_elip] = elip_curve(a, b0)
syms t
% parameterization of x(t), y(t) -- from Horn pg. 4
b = b0;
x_elip = a*cos(t);
y_elip = b*sin(t);

eccentricity = sqrt(1 - (a/b)^2);

% k_elip = (a*b) / (a^2*(sin(t))^2 + b^2*(cos(t))^2)^(3/2);

% expression to be integrated -- Horn pg. 5
expr = (1 / (1 - eccentricity^2*(sin(t))^2)^(5/2));
% expr_2 = k_elip^2 * sqrt(diff(x_elip,t,1)^2 + diff(y_elip,t,1)^2);

E_elip = (a^2/b^3) * vpa(int(expr, 0, pi/2));
end
