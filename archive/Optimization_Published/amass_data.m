function amass_data( bases, curves )
%AMASS_DATA Ammass the data resulting from the 
%           optimizations.

global divs;

global perc_diff;
global base_curvatures;
global new_curvatures;
global diff_curvatures;
global pdiff_curvatures;
global base_derivative1;
global base_derivative2;
global new_curve_derivative1;
global new_curve_derivative2;

[m n] = size(curves);

perc_diff = [];
base_curvatures = [];
new_curvatures = [];
diff_curvatures = [];
pdiff_curvatures = [];
base_derivative1 = [];
base_derivative2 = [];
new_curve_derivative1 = [];
new_curve_derivative2 = [];

for i=1:m
    b = bases(i,1:divs);
    p = curves(i,1:divs);

    % Compute the percent change of the radius's
    perc_diff = cat(1, perc_diff, 100*(p-b)./b);
    
    % Compute the derivatives of the old curve
    b1 = deriv(b, 2,divs-1);
    b2 = deriv(b1,3,divs-1);
    
    % Compute the derivatives of the new curve
    p1 = deriv(p, 2,divs-1);
    p2 = deriv(p1,3,divs-2);
    
    curve = 3:divs-2;
    curve0 = 3:divs-2;
    dcurve = 3:divs-2;
    pdcurve = 3:divs-2;
    
    % Sum the difference in the square of the curvatures
    for j=3:divs-2
    
        k  = kappa(p(j),p1(j),p2(j));
        k0 = kappa(b(j),b1(j),b2(j));
        
        curve(j) = k;
        curve0(j) = k0;
        dcurve(j) = k-k0;
        pdcurve(j) = 100*(k-k0)/k0;
    end
    
    base_curvatures = cat(1, base_curvatures,curve0);
    new_curvatures = cat(1, new_curvatures,curve);
    diff_curvatures = cat(1, diff_curvatures, dcurve);
    pdiff_curvatures = cat(1, pdiff_curvatures, pdcurve);
    base_derivative1 = cat(1,base_derivative1,b1);
    base_derivative2 = cat(1,base_derivative2,b2);
    new_curve_derivative1 = cat(1,new_curve_derivative1,p1);
    new_curve_derivative2 = cat(1,new_curve_derivative2,p2);

end

end

