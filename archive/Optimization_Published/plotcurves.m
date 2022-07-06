function plotcurves( num, final )
%PLOTCURVES Summary of this function goes here
%   Detailed explanation goes here

global theta;
global base;

% Use the minor axis as a number for the figure
figure(num);

% Plot both the final in blue and base in red
polar(theta,final);
hold on
polar(theta,base,'r');
hold off

end

