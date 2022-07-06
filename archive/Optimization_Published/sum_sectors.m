function [ area ] = sum_sectors( radius )
%SUM_SECTORS Ingtegral approximation by summing the
%       areas of the sectors determined from an array
%       of radius's and theta's.

global dtheta;
global divs;

% Initialize the area at zero
area = 0;

% Calculate the area of the ellipse by summing area of sectors
for i=3:divs-2
    r = radius(i);
    sector = (r.^2)*cos(dtheta/2)*sin(dtheta/2); 
    area = area + sector; 
end
