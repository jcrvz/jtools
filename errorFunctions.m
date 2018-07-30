function [SSE,R2,MSE,RMSE] = errorFunctions(xMeasured,yMeasured,modelToFit)
    % modelToFit must have the form: y <- @(p,x) where p is the vector of
    % unknown parameters and x is the vector independent variables (known)
    %
    % This function returns following error functions:
    % R2    - the R-squared function (R^2)
    % SSE   - the sum of squares due to error function
    % MSE   - the mean squared error function
    % RMSE  - the root mean squared error function
    
    % Determine the number of data
    N = length(xMeasured);
    
    % Find the total sum of squares (TSS)
    TSS = sum((yMeasured - mean(yMeasured)).^2);
    
    % Create the sum of squares due to error (SSE) function
    SSE = @(p) sum((yMeasured - modelToFit(p,xMeasured)).^2);
    
    % Create the R-squared function (R^2)
    R2 = @(p) 1 - SSE(p)/TSS;
    
    % Create the mean squared error (MSE) function
    MSE = @(p) SSE(p)/N;
    
    % Create the root mean squared error (RMSE) function
    RMSE = @(p) sqrt(MSE(p));
end