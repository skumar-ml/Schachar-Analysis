function [x_param, y_param] = fourier_fit(X_data, Y_data)
%FOURIER_CURVE 
% Given an age, returns the parameterized syms equation for 
% the curve (in cartesian coordinates)
    

    theta = atan2(Y_data, X_data);
    rho = sqrt(X_data.^2 + Y_data.^2);
    
%     figure; hold on;
%     scatter(X_data, Y_data)
%     scatter(theta, rho)

    fun = @(x, t)x(1)*cos(0) + x(2)*cos(t) + x(3)*cos(2*t) + x(4)*cos(3*t) + x(5)*cos(4*t) + x(6)*cos(5*t) + x(7)*cos(6*t) + x(8)*cos(7*t) + x(9)*cos(8*t) + x(10)*cos(9*t) + x(11)*cos(10*t);
    x0 = ones(1,11);
    
%     x = lsqcurvefit(fun, x0, acos(X_data ./ max(X_data)), Y_data);
    coef = lsqcurvefit(fun, x0, theta, rho);
    
    syms t
    
    sym_eq = sym(fun(coef, t));
    
    x_param = sym_eq * cos(t);
    y_param = sym_eq * sin(t);
    
%     fp = fplot(x_param, y_param, [pi/2, 3*pi/2], 'LineWidth', 3);
%     X_fourier = fp.XData; Y_fourier = fp.YData;
%     scatter(atan2(Y_fourier, X_fourier), sqrt(X_fourier.^2 + Y_fourier.^2))
%     legend("Raw data (Cartesian)", "Raw data (Polar)", "Fourier (Cartesian)", "Fourier (Polar)")
%     figure

end