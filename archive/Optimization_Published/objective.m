function [ value ] = objective( p )
%OBJECTIVE function for constrained optimization of ellipse
%          We want to minimize the difference between the
%          curvatures.

global base;
global base1;
global base2;
global divs;
global dtheta;

% Store the derivatives of the baseline curve in local variables
b  = base;
b1 = base1;
b2 = base2;

% Compute the derivatives of the new curve
p1 = deriv(p, 2,divs-1);
p2 = deriv(p1,3,divs-2);

% Iinitialize the result of the objective
value = 0;

% Sum the difference in the square of the curvatures
for i=3:divs-2
    
    k  = kappa(p(i),p1(i),p2(i));
    k0 = kappa(b(i),b1(i),b2(i));
    term1 = (k - k0)^2;
    term2 = ((p(i))^2 + (p1(i))^2)^(1/2);
    product = term1*term2;
    sector = (product.^2)*cos(dtheta/2)*sin(dtheta/2); 
    value = value + sector;

end


end