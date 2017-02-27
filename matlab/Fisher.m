% FIM assuming unweighted least squares with x0 at end of params 

function FIM = Fisher(tspan,params,xfcn,yfcn,data_case)
% The plan--
%  recall the measurement y
%  for each time point in t we need a sensitivity vector of dy_i/dp's
%   - note we'll use 1% perturbations in the parameter to calculate these
%  the design matrix X is formed by taking a column vector of the transpose
%  of these vectors. (See Landaw notes II-10, Biomath 270)
%  FIM = X'X (See Landaw II-12)
%
% Note: we won't actually calculate the sensitivity vectors for each time 
% point one at a time; instead we'll calculate the sensitivity for each 
% parameter, for all time pts simultaneously, then fill in the columns for X
% one by one.

% design matrix
% X = zeros(length(tspan)-1,length(params));  %we do tspan-1 if the initial conditions 
% (for this version) are fixed, so the sensitivity at t_0 will be 0.
X = [];
for j=[1,2,6,7,8,9] %instead of 1:length(params), because we fixed rest of the parameters
    params1 = params;
    params2 = params;
    params1(j) = 1.01*params(j);
    params2(j) = 0.99*params(j);
    x0fcn = @(data,params) [1-(3*data(1)/params(end)); 2*data(1)/params(end); ...
    data(1)/params(end); 1; 1/params(8); 0;0];
    x01 = x0fcn(data_case,params1);
    x02 = x0fcn(data_case,params2);
    [t x1] = ode45(xfcn,tspan,x01,[],params1);
    [t x2] = ode45(xfcn,tspan,x02,[],params2);
    X = [X, (yfcn(x1,params1) - yfcn(x2,params2))./(0.02*params(j))];
    %this fills in the jth column of the design matrix with the sensitivities to parameter j 
    %at each time point.
end

% COV from FIM
FIM = X'*X;

%=== end ===%