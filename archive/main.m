clear; close all; clc;
load strenkA_anterior.mat;

ratio = 0.05 : 0.05 : 1;
figure

hold on
plot(ratio, E_ron_total)
plot(ratio, E_elip_total)
legend("Chien", "Ellipse")
title("Energy, Strenk A Anterior")
xlabel("Ratio (b0/a)")