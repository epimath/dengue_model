% This function calculates the RSS for the Dengue Model
function RSS = dengueRSS(DataTimes,DataVals,yfcn,x0fcn,params,paramsfixedfcn)
params = abs(params);

% Fixed Params
params = paramsfixedfcn(params);

% Initial Conditions
x0 = x0fcn(DataVals,params);
% x0 = zeros(1,7);
% x0(3) = DataVals(1)/params(end); %0.000001;
% x0(2) = 10*x0(3);
% x0(1) = 1 - x0(3) - x0(2);
% %x0(4) = 2*(params(6) - params(7)*params(8))/params(6);
% x0(4) = 1; 
% x0(5) = x0(4)/params(8); 


% Simulate the ODE
[t,x] = ode45(@denguemodel,DataTimes,x0,[],params);

% Calculate the output/measurement equation
y = yfcn(x,params);  

% Calculate the weights
W = diag(ones(1,length(y))./ mean2(y));

% Calculate the weighted sum of squared residuals (least squares)
RSS = ((DataVals' - y))'*W*((DataVals' - y));


%=== end ===%