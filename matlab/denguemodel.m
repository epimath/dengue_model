% ODE function for the dengue model
function dxdt = denguemodel(t,x,param)
dxdt = zeros(7,1);

%Sh Eh Ih A Sm Em Im
%human units: fractions
%larvae units: fractions of carrying capacity
%adult mosquito units: scaled to carrying capacity & maturation rate

betamh = param(1); %transmission rate from mosquito to human
betahm = param(2); %transmission rate from human to mosquito
alpha = param(3);  %intrinsic incubation rate
eta = param(4);    %recovery rate
gamma = param(5);  %extrinsic incubation rate 
k = param(6);      %case reporting rate
mua = param(7);    %aquatic death rate
mum = param(8);    %mosquito death rate
%mu = param(9);    %human death rate

% Note: we've dropped mu for humans for now---the outbreak is half a year, 
% so we can reasonably neglect birth-death dynamics
dxdt(1) = - betamh * x(1) * x(7);
dxdt(2) = betamh * x(1) * x(7) -alpha * x(2);
dxdt(3) = alpha * x(2) - eta * x(3);
dxdt(4) = k*(x(5) + x(6) + x(7)) * (1 - x(4)) - mua * x(4);
dxdt(5) = x(4) - betahm * x(5) * x(3) - mum * x(5);
dxdt(6) = betahm * x(5) * x(3) - gamma * x(6) - mum * x(6);
dxdt(7) = gamma * x(6) - mum * x(7);

%SIR Only
% dxdt(1) = - betamh * x(1) * x(3);
% dxdt(2) = betamh * x(1) * x(3) -alpha * x(2);
% dxdt(3) = alpha * x(2) - eta * x(3);

%=== end ===%