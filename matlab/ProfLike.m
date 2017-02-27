% Profile Likelihood Generator

function profile = ProfLike(params,profparam,fitter)
% Definitions
%   params = point in parameter space from which to profile (i.e. the parameter estimates)
%   profparam = index of the parameter to be profiled
%   fitter = this is a customized fminsearch that takes two arguments:
%     params and paramsfixedfcn, which tell it the starting parameters and
%     fixes the profiled parameter. Everythng else (data, ICs, likelihood
%     function, etc.) is fixed for the entire profile so is set outside when 
%     the fitter is defined.
%     e.g. fitter = @(params,paramfixedfcn) fminsearch(@(p) siwrML(times,p,paramfixedfcn,data,x0fcn,yfcn),params,optimset('MaxFunEvals',5000,'MaxIter',5000))


% Setup
per = 0.1; %0.5
numpoints = 8; %10

% Profile
profrangeDown = linspace(params(profparam), params(profparam)*per,numpoints)'; 
profrangeUp = linspace(params(profparam), params(profparam)/per,numpoints)';
% split into up and down so we can use last fitted value as starting value for next run
profrange = [profrangeDown profrangeUp];
currfvals = [];
currparams = [];
currflags = [];
for i=1:2
    paramstemp = params;
    for j = 1:numpoints
        [i j] %track progress
        if profparam==length(params) % got to be a nicer way to do this without an if statement, but think of it later
            paramfixedfcn = @(p) [p(1:profparam-1) profrange(j,i)]; % for the last parameter in the array
        else
            paramfixedfcn = @(p) [p(1:profparam-1) profrange(j,i) p(profparam+1:end)];
        end
        [paramstemp, fvaltemp, flagtemp] = fitter(paramstemp,paramfixedfcn);
        paramstemp = paramfixedfcn(abs(paramstemp));
        currfvals = [currfvals; fvaltemp];
        currflags = [currflags; flagtemp];
        currparams = [currparams; paramstemp];
    end
end
profile = [flipud([profrangeDown currfvals(1:numpoints) currflags(1:numpoints) currparams(1:numpoints,:)]);...
   [profrangeUp currfvals(numpoints+1:end) currflags(numpoints+1:end) currparams(numpoints+1:end,:)]];

% Old version from when it wasn't a function
% profiles = []; %SLAFD U+1F47B
% for i=1:2%length(params)
%     profrange = [linspace(params(i)/factor, params(i),numpoints/2)'; linspace(params(i), params(i)*factor,numpoints/2)'];
%     currfvals = [];
%     currparams = [];
%     currflags = [];
%     for j = 1:length(profrange)
%         [i j] %track progress
%         if i==length(params) % got to be a nicer way to do this without an if statement, but think of it later
%             paramfixedfcn = @(p) [p(1:i-1); profrange(j)];
%         else
%             paramfixedfcn = @(p) [p(1:i-1); profrange(j); p(i+1:end)];
%         end
%         [paramstemp, fvaltemp, flagtemp] = fminsearch(@(p) siwrML(times,p,paramfixedfcn,data,x0fcn,yfcn),params,optimset('MaxFunEvals',5000,'MaxIter',5000));
%         paramstemp = abs(paramstemp);
%         currfvals = [currfvals; fvaltemp];
%         currflags = [currflags; flagtemp];
%         currparams = [currparams; paramstemp'];
%     end
%     profiles(:,:,i) = [profrange currfvals currflags currparams];
%     
%     figure(i)
%         hold on
%         plot(profrange,currfvals,'k.-','LineWidth',2)
%         plot(params(i),fval,'ro','LineWidth',2)


%=== end ===%


