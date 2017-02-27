%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Main code:        %
% run and fit dengue model %
% Author: Marisa Eisenberg %
% Email: marisae@umich.edu %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Note:
% denguemodel.m (model equations), dengueRSS.m (cost fucntion), 
% Fisher.m (Fisher information matrix), ProfLike.m (profile likelihood) are called in this script


%% start
clear

%% Data

% Dengue cases each week
data_case=[6.0, 4.0, 25.0,20.0, 29.0, 37.0, 71.0, 54.0, 77.0, 52.0, 70.0, 95.0, 70.0, 66.0,78.0, 53.0, 67.0, 52.0, 62.0, 33.0, 13.0, 14.0, 5.0, 0.0, 2.0, 1.0,1.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0];

% Data Times are in Days
data_times = 0:7:7*length(data_case)-1;


%% Model Setup: Initial Conditions, Parameters, Measurement Equations

% Parameter guess (already set to the fitted values here)
betamh =14.1530277;
betahm =0.0300024024; 
alpha = 0.14; 
eta = 0.2;
gamma = 0.1; 
k = 2.02896263;
mua = 4.18336111;
mum = 0.316770764;
reph = 1546.73938; 
params = [betamh betahm alpha eta gamma k mua mum reph];


% Fixed Parameters - made this a function so that if we change the setup
% later it will automatically pass to all the other spots where we set
% fixed params. If you want to fit all parameters, just make this return params.
paramsfixedfcn = @(params) [params(1:2) alpha eta gamma params(6:end)];


% Initial Conditions - made this a function so that if we change the setup
% later it will automatically pass to all the other spots where we set
% initial conditions. 
% Note: Eh = 2*Ih, A = CC = 1, Sm = SS with A, Em = Im = 0
x0fcn = @(data,params) [1-(3*data(1)/params(end)); 2*data(1)/params(end); ...
    data(1)/params(end); 1; 1/params(8); 0;0];


% Measurement Equations
yfcn = @(x,params) 7*(cumtrapz(params(end)*params(3)*x(:,2))+6/7) - 7*[0; cumtrapz(params(end)*params(3)*x(1:end-1,2))+6/7];        %weekly incidence (i.e. cumulative for the week)
% Note: assumes weekly time points!!  fix this later to make flexible
% also note: 7* is to make it spacings of 7 instead of 1 between pts (note linear so this is legal)
% the +6/7 is to start the integral with 6 cumulative cases so that we match the initial data (with correction for the *7). 


%% Parameter Estimation
ParamEsts = fminsearch(@(params) dengueRSS(data_times,data_case,yfcn,x0fcn,params,paramsfixedfcn),params,optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000));
ParamEsts = abs(ParamEsts);
ParamEsts = paramsfixedfcn(ParamEsts);
%ParamEsts = [betamh betahm alpha eta gamma k mua mum reph];
%params = ParamEsts;


%% Simulate Model
% Initial Conditions
x0 = x0fcn(data_case,ParamEsts);

% Simulate the model
[t_est,x_est] = ode45(@denguemodel,data_times,x0,[],ParamEsts);

% Calculate y (measurement)
y_est = yfcn(x_est,ParamEsts);


%% Plot data and fitted curve

figure(1)
    set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
    hold on
    plot(t_est,y_est,'k','LineWidth',2);
    plot(data_times,data_case,'o');
    ylabel('Dengue Cases');  
    xlabel('Time (days)');
    
figure(2)
    set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
    hold on
    plot(t_est,x_est(:,4:end),'LineWidth',2);
    legend('Baby mosquitoes','Susceptible adults','Exposed adults','Infected adults');  
    xlabel('Time (days)');


%% Calculate & Print R0

Sm_DFE = (ParamEsts(6) - ParamEsts(7)*ParamEsts(8))/(ParamEsts(6)*ParamEsts(8))
Mu = 0; %placeholder in case we add it back in later

R0 = (sqrt(Sm_DFE*ParamEsts(3)*ParamEsts(2)*ParamEsts(1))*ParamEsts(5))...
    /sqrt((ParamEsts(3)+Mu)*(ParamEsts(4)+Mu)*ParamEsts(8)*ParamEsts(5)*(ParamEsts(5)+ParamEsts(8)));

%% Calculate Fisher information matrix
fisherm = Fisher(data_times,paramEsts,@denguemodel,yfcn,data_case);
f_rank = rank (fisherm);

   
%% Calculate Profiles
times = t_est;
data = y_est;
profiles = []; %SLAFD U+1F47B
fitter = @(params,paramfixedfcn) fminsearch(@(p) dengueRSS(data_times,data_case,yfcn,x0fcn,p,paramsfixedfcn),paramEsts,optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000));
%threshold = (fval/(length(data) - length(paramests)))*chi2inv(0.95,length(paramests)); %need to add fval if it's fitted to actual data (or noisy data); 

% Plot profiles
for i=1:length(paramEsts)
    profiles(:,:,i) = ProfLike(paramEsts,i,fitter);
end

paramlist = {'betamh', 'betahm', 'alpha', 'eta', 'gamma', 'k', 'mua', 'mum', 'reph'};
for i=1:length(params)
    %plot profile
    figure(i)
    set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
    hold on
    plot(profiles(:,1,i),profiles(:,2,i),'k','LineWidth',2)
    %plot(paramEsts(i),fval,'r*','LineWidth',2)
    %plot(profiles(:,1,i),threshold*ones(size(profiles(:,1,i))),'k--')
    xlabel(paramlist{i})
    ylabel('Cost Function Value')
    
    %plot parameter relationships
    figure(10+i)
    set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
    hold on
    plot(profiles(:,1,i),profiles(:,4:end,i),'LineWidth',2)
    plot(params(i),params,'r*')
    xlabel(paramlist{i})
    ylabel('Estimated Parameter Value')
    legend(paramlist)
end

%=== end ===%


