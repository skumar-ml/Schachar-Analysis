clear; close all; clc;
%% Input data
data_file = "Human_20 yo RIEB15-1632_OD_data.xls";
age = 20;

%% Read & Process Data

% Read data
data_path = strcat("data/", data_file);

% Reads all X/Y data from correct sheet and stores in matrix M. 
% First column = x, Second column = y

M = readmatrix(data_path, 'Sheet', 'Centered and Aligned', 'Range', 'A:B');
figure; scatter(M(:,1), M(:,2)); title("Original plot of data"); xlim([-5,5]); ylim([-5,5]); % Plot to confirm

% Bottom in original data, is indented from suture, so I
% replicate from top side to smoothen out data and make consistent

post = M(M(:, 1) > 0, :); % Filter out anterior
figure; scatter(post(:,1), post(:,2)); title("Posterior Raw"); % Plot to check

post_top = post(post(:,2) > 0, :); % Filter out top of posterior
post_bot = [post_top(:,1), -1*post_top(:,2)]; % Flip across x-axis

post_new = cat(1, post_top, post_bot); % Concat to form new anterior
figure; scatter(post_new(:,1), post_new(:,2)); title("Posterior w/ Fixed Suture") % Plot to check

X = post_new(:,1); Y = post_new(:,2); % Draw out X, Y from anterior data

% Process data for algorithm (places anterior on top, optic axis is x-axis)
X_data = Y; Y_data = X;
figure; scatter(X_data, Y_data, 6); hold on;

%% Fourier 

% Find parameterized x, y of Fourier Curve
% [x_fourier, y_fourier] = fourier_curve(age);

data_top = M(M(:,2) > 0, :);
data_bottom = [data_top(:,1), -1*data_top(:,2)];

fourier_data = cat(1, data_top, data_bottom);
[x_fourier, y_fourier] = fourier_fit(fourier_data(:,1), fourier_data(:,2));

temp = x_fourier; x_fourier = y_fourier; y_fourier = temp;

% Obtain x,y coords for posterior
fourier_bounds = [-pi/2, pi/2];
fp = fplot(x_fourier, y_fourier, fourier_bounds); X_fourier = fp.XData; Y_fourier = fp.YData;
%plot(X_fourier, Y_fourier, 'LineWidth', 2);

a_post_fourier = max(X_fourier);

%% Chien - fit to raw
syms t

% Get anterior a

ant = M(M(:, 1) < 0, :); % Filter out anterior

ant_top = ant(ant(:,2) > 0, :); % Filter out top of anterior
ant_bot = [ant_top(:,1), -1*ant_top(:,2)]; % Flip across x-axis

ant_new = cat(1, ant_top, ant_bot); % Concat to form new anterior

X_anterior = ant_new(:,1); Y_anterior = ant_new(:,2); % Draw out X, Y from anterior data

% Process data for algorithm (places anterior on top, optic axis is x-axis)
X_ant = Y_anterior; Y_ant = -X_anterior;

% Chien analysis
b0 = max(Y_data);
a = max(X_data) + 0.0001;
a_ant = max(X_ant)+ 0.0001; % add epsilon for numerical stability

[b1, b3] = findChienCoefficients(X_data', Y_data', a, b0);

x_chien = a_ant*sin(t);
y_chien = (b0 + b1*t^2 + b3*t^4)*cos(t);

chien_bounds = [-pi/2, pi/2];
fp = fplot(x_chien, y_chien, chien_bounds); X_chien = fp.XData; Y_chien = fp.YData;
%plot(X_chien, Y_chien, 'LineWidth', 2);


%% Forbes
% format data to forbes specs
Y_forbes = -1*Y_data + max(Y_data); %figure; scatter(X_data, Y_forbes);

% fit forbes to data
syms rho;
[forbes_eq, Y_forbes_raw, A] = forbes(X_data', Y_forbes', 6);
forbes_reformat = -1*forbes_eq + eval(subs(forbes_eq, rho, a));
forbes_eq = forbes_reformat;

fp = fplot(rho, forbes_eq, [min(X_data), max(X_data)]); X_forbes = fp.XData; Y_forbes = fp.YData;
%Y_forbes = -Y_forbes + max(Y_forbes); % Reconvert to standard format (matching other graphs)
%plot(X_forbes, Y_forbes, 'LineWidth', 2);

forbes_eq = subs(forbes_eq, rho, t);
% Note -- the t in the forbes equation is cartesian! Stands for x (not
% theta)
%% Ellipse
x_elip = a*cos(t); % in mm
y_elip = b0*sin(t); % in mm

elip_bounds = [0 pi];
fp = fplot(x_elip, y_elip, elip_bounds); X_elip = fp.XData; Y_elip = fp.YData;
%plot(X_elip, Y_elip, 'LineWidth', 2);

legend("Raw", "Fourier", "Chien", "Forbes", "Ellipse")

%% Metrics
zone = 3; % optical zone [-zone, +zone]
offset_chien = abs(chien_bounds(1)) - asin(zone/a); % polar - difference to come in from edges
offset_elip = -1* (abs(elip_bounds(1)) - acos(zone/a));
offset_fourier = abs(fourier_bounds(1)) - asin(zone/a_post_fourier)

% Find curvature
k_chien = findCurvature(x_chien, y_chien, chien_bounds(1)+offset_chien, chien_bounds(2)-offset_chien);
k_elip = findCurvature(x_elip, y_elip, elip_bounds(1)+offset_elip, elip_bounds(2)-offset_elip);
k_fourier = findCurvature(x_fourier, y_fourier, fourier_bounds(1)+offset_fourier, fourier_bounds(2)-offset_fourier);
k_forbes = findCurvature(t, forbes_eq, -zone, zone)

% Plot curvature
figure; hold on;
fplot(a_ant*sin(t), abs(k_chien), chien_bounds);
fplot(a_ant*sin(t-pi/2), abs(k_elip), elip_bounds);
fplot(a_post_fourier*sin(t), abs(k_fourier), fourier_bounds);
fplot(t, abs(k_forbes));
legend("Chien", "Ellipse", "Fourier", "Forbes"); title("Curvature")

% Find smoothing energy (integral of derivative of curvature squared)
smth_chien = eval(vpaintegral(diff(k_chien, t, 1) ^ 2, chien_bounds(1)+offset_chien, chien_bounds(2)-offset_chien));
smth_elip = eval(vpaintegral(diff(k_elip, t, 1) ^ 2, elip_bounds(1)+offset_elip, elip_bounds(2)-offset_elip));
smth_fourier = eval(vpaintegral(diff(k_fourier, t, 1) ^ 2, fourier_bounds(1)+offset_fourier, fourier_bounds(2)-offset_fourier));
smth_forbes = eval(vpaintegral(diff(k_forbes, t, 1) ^ 2, -zone, zone))

% Find bending energy
[bendE_chien, firstD_chien, expr_chien] = findBendingEnergy(x_chien, y_chien, chien_bounds(1)+offset_chien, chien_bounds(2)-offset_chien);
[bendE_elip, firstD_elipAnt, expr_elip] = findBendingEnergy(x_elip, y_elip, elip_bounds(1)+offset_elip, elip_bounds(2)-offset_elip);
[bendE_fourier, firstD_fourierAnt, expr_fourier] = findBendingEnergy(x_fourier, y_fourier, fourier_bounds(1)+offset_fourier, fourier_bounds(2)-offset_fourier);
[bendE_forbes, firstD_forbes, expr_forbes] = findBendingEnergy(t, forbes_eq, -zone, zone)

figure
% Mean/Variance of RoC - variance found numerically
meanROC_chien = abs(1/((chien_bounds(2)-offset_chien) - (chien_bounds(1)+offset_chien)) * eval(vpaintegral(k_chien, chien_bounds(1)+offset_chien, chien_bounds(2)-offset_chien)));
fp = fplot(abs(k_chien), [chien_bounds(1)+offset_chien, chien_bounds(2)-offset_chien], 'MeshDensity', 200);
yROC_chien = fp.YData;
varROC_chien = var(yROC_chien);
valROC_chien = abs(eval(subs(1/k_chien, t, chien_bounds(1)+offset_chien)))


meanROC_elip = abs(1/((elip_bounds(2)-offset_elip) - (elip_bounds(1)+offset_elip)) * eval(vpaintegral(k_elip, elip_bounds(1)+offset_elip, elip_bounds(2)-offset_elip)));
fp = fplot(abs(k_elip), [elip_bounds(1)+offset_elip, elip_bounds(2)-offset_elip], 'MeshDensity', 200);
yROC_elip = fp.YData;
varROC_elip = var(yROC_elip);
valROC_elip = abs(eval(subs(1/k_elip, t, elip_bounds(1)+offset_elip)))

meanROC_fourier = abs(1/((fourier_bounds(2)-offset_fourier) - (fourier_bounds(1)+offset_fourier)) * eval(vpaintegral(k_fourier, fourier_bounds(1)+offset_fourier, fourier_bounds(2)-offset_fourier)));
fp = fplot(abs(k_fourier), [fourier_bounds(1)+offset_fourier, fourier_bounds(2)-offset_fourier], 'MeshDensity', 200);
yROC_fourier = fp.YData;
varROC_fourier = var(yROC_fourier);
valROC_fourier = abs(eval(subs(1/k_fourier, t, fourier_bounds(1)+offset_fourier)))

meanROC_forbes = abs(1/(2*zone) * eval(vpaintegral(k_forbes, -zone, zone)));
fp = fplot(abs(k_forbes), [-zone, zone], 'MeshDensity', 200);
yROC_forbes = fp.YData;
varROC_forbes = var(yROC_forbes);
valROC_forbes = abs(eval(subs(1/k_forbes, t, -zone)))

% Fit

data = [X_data, Y_data];
data_fit = data(-3 < data(:,1), :);
data_fit = data_fit(data_fit(:,1) < 3, :);
X_fit = data_fit(:, 1); Y_fit = data_fit(:, 2);
figure; scatter(X_fit, Y_fit)

fit_forbes = getFit(X_fit, Y_fit, forbes_eq);
fit_elip = getFit(acos(X_fit ./ a), Y_fit, y_elip);
fit_chien = getFit(asin(X_fit ./ a), Y_fit, y_chien);
fit_fourier = getFit(asin(X_fit ./ a_post_fourier), Y_fit, y_fourier);


