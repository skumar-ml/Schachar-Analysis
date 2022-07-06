clear; close all;

%% Fourier Curve
syms t
age = 20;
% Find parameterized x, y of Fourier Curve
[x_fourier, y_fourier] = fourier_curve(age);

% Obtain x,y coords
[Y_fourier, X_fourier] = fplot(x_fourier, y_fourier, [0, 2*pi]);
plot(X_fourier, Y_fourier)

%% Fitting Chien curve - Posterior
fourier_bounds = [-pi/2, pi/2]
[Y_fourierPost, X_fourierPost] = fplot(x_fourier, y_fourier, fourier_bounds);
figure
plot(X_fourierPost, Y_fourierPost)

% Posterior a/b
b0_post = max(Y_fourierPost);
a_post = max(X_fourierPost);

[b1_post, b3_post] = findChienCoefficients(X_fourierPost, Y_fourierPost, a_post, b0_post);

x_chienPost = a_post*sin(t);
y_chienPost = (b0_post + b1_post*t^2 + b3_post*t^4)*cos(t);

chien_bounds = [-pi/2, pi/2]
[X_chienPost, Y_chienPost] = fplot(x_chienPost, y_chienPost, chien_bounds);

hold on
plot(X_chienPost, Y_chienPost)


%% Ellipse - Posterior
x_elipPost = a_post*cos(t);
y_elipPost = b0_post*sin(t);

elip_bounds = [0 pi];
[X_elipPost, Y_elipPost] = fplot(x_elipPost, y_elipPost, elip_bounds);
plot(X_elipPost, Y_elipPost)

legend('Fourier', 'Chien', 'Ellipse')
title('Posterior Curve')

%% Posterior Energy (Fourier, Chien, Ellipse)
[E_fourierPost, k_fourierPost, var_fourierPost]  = findEnergy(x_fourier, y_fourier, fourier_bounds(1), fourier_bounds(2))

[E_chienPost, k_chienPost, var_chienPost] = findEnergy(x_chienPost, y_chienPost, chien_bounds(1), chien_bounds(2))

[E_elipPost, k_elipPost, var_elipPost] = findEnergy(x_elipPost, y_elipPost, elip_bounds(1), elip_bounds(2))

% Graph of curvatures
figure
fplot(t+pi/2, abs(k_fourierPost), fourier_bounds)
hold on
fplot(t+pi/2, abs(k_chienPost), chien_bounds)
fplot(t, abs(k_elipPost), elip_bounds)

legend("Fourier", "Chien", "Ellipse")
title("Curvature for Posterior")

