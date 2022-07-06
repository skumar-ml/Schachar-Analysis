function plots(base,curves)
%PLOTS Summary of this function goes here
%   Detailed explanation goes here

global perc_diff;
global pdiff_curvatures;

% Set up figures of results
aspectratio = zeros(8,1);
for i=3:10
    aspectratio(i-2,1) = i/10;
end

y1 = pdiff_curvatures(1:8,103);
y2 = perc_diff(1:8,103);

YMatrix = cat(2,y1,y2);
createfigure1(aspectratio, YMatrix);

figure(2);
polar(theta,curves(8,1:105));
hold on;
polar(theta,base(8,1:105),'r');
for i=1:7
    polar(theta,curves_secant((7-i),1:105));
    polar(theta,base_curves((7-i),1:105),'r');
end
hold off

end