function ellipse_setup ( a, b, divisions)
%RHO_0 Output the value of the radius of the ellipse
%   with semi-axis's a and b and degree theta

clearvars -global;
global base;
global theta;
global divs;
global dtheta;
global barclength;
global initial;
global ai;
global bi;

global base1;
global base2;

% Add divisions for calculating derivatives at boundaries
divs = divisions+5;

% Caclulate the angle of each division
i = pi/2/divisions;
dtheta = i;


% Calculate the baseline radius
for j=1:divs
    term1 = (a^2)*((sin((j-3)*i))^2);
    term2 = (b^2)*((cos((j-3)*i))^2);
    denom = sqrt(term1+term2);
    numer = a*b;
    base(j) = numer./denom;
    theta(j) = (j-3)*i;
end

% Calculate the initial guess radius's
ai = a*(1.02);
bi = b*(0.98);
for j=1:divs
    term1 = (ai.^2)*((sin((j-3)*i)).^2);
    term2 = (bi.^2)*((cos((j-3)*i)).^2);
    denom = sqrt(term1+term2);
    numer = ai*bi;
    initial(j) = numer./denom;
end

% Calculate the derivatives of the baseline
base1 = deriv(base ,2,divs-1);
base2 = deriv(base1,3,divs-2);

% Initialize base arc length to zero
barclength = 0;

% Calculate the baseline arc length
for j=3:divs-3
    x1 = base(j)*cos(theta(j));
    y1 = base(j)*sin(theta(j));
    x2 = base(j+1)*cos(theta(j+1));
    y2 = base(j+1)*sin(theta(j+1));
    l = ((x2-x1)^2 + (y2-y1)^2)^(1/2);
    barclength = barclength + l;
end

% % Graph and display the baseline and derivatives
% figure(11);
% polar(theta,base);

% % Plot the baseline and derivatives
% figure(22);
% t = 1:divs;
% subplot(3,1,1);
% plot(t,base);
% subplot(3,1,2);
% plot(t,base1,'r');
% subplot(3,1,3);
% plot(t,base2,'g');

end