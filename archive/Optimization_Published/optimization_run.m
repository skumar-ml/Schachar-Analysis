function new_curves = optimization_run()
% Runs the optimization, iterating over aspect ratios

% The resulting final curve after minimization
new_curves = [];

% Start with the default options
options = optimset;
% Modify options setting for the problem.  
%options = optimset(options,'Display', 'iter');
options = optimset(options,'MaxFunEvals', 500000);
options = optimset(options,'MaxIter', 10000);
options = optimset(options,'Algorithm', 'interior-point');
%options = optimset(options,'Diagnostics', 'on');

% Iterate from an aspect ratio of .3 to 1
for b=3:10
    
    % Indicate which run we are starting
    str = sprintf('Starting the run on ellipse with aspect ratio %d', (b/10));
    disp(str);
    
    % Set up the mimization of 10 by b ellipse with 100 divisions
    ellipse_setup(10,b,100);
    global initial;
    
    % Call the optimization tool, save new curve in x
    [x,fval,exitflag,output,lambda,grad,hessian] = ...
    fmincon(@objective,initial,[],[],[],[],[],[],@constraints,options);
    
    % Plot the result of base ellipse vs. new curve
    plotcurves(b,x);
    disp(x);
    
    % Store the new curve result by adding to the new curve array
    new_curves = cat(1,new_curves,x);
    
    % Indicate which run we just finished
    str = sprintf('Done with run on ellipse with aspect ratio %d', (b/10));
    disp(str);
end

end
