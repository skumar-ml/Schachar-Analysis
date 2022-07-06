%% Horn's method - reproduced
syms y
m = (1/sqrt(2))^(2); % Matlab uses m = k^2 (Horn uses k for his integral)
c = 0.8346; % the constant value Horn uses for his example case

% Horn provides an equation relating y to theta
y_min = 0;  % theta = pi/2
y_max = 2*c; % theta = 0

x = sqrt(2)*c*(2*ellipticE(acos(y/(2*c)),m) - ellipticF(acos(y/(2*c)),m)); % Horn's cartesian equation
x_inv = -sqrt(2)*c*(2*ellipticE(acos(y/(2*c)),m) - ellipticF(acos(y/(2*c)),m));

[X1, Y1] = fplot(x);
[X2, Y2] = fplot(x_inv);
figure
plot(Y1, X1) % to switch major and minor axis, just swap X and Y
hold on
plot(Y2, X2)
title("Horn's method, reproduced")

figure
[X_horn, Y_horn] = fplot(x, [y_min y_max]);
plot(X_horn, Y_horn)
title("Horn's method, entire range of theta (from pi/2 to 0)")