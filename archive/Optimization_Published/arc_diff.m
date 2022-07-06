function percdiff = arc_diff( r )
%ARC_LENGTH Calculate the percent difference in arc length
%           between the given set of radius's r and the baseline.

global barclength;
global divs;
global theta;

% Set length initially to zero
l = 0;

% Calculate the baseline arc length
for j=3:divs-3
    x1 = r(j)*cos(theta(j));
    y1 = r(j)*sin(theta(j));
    x2 = r(j+1)*cos(theta(j+1));
    y2 = r(j+1)*sin(theta(j+1));
    section = ((x2-x1)^2 + (y2-y1)^2)^(1/2);
    l = l + section;
end

% Calculate the percent difference from 
percdiff = (l - barclength)/(barclength);

end

