function [ID, IDn] = ID_name(name,theta,nx,ny,nvx,nvy,RK_stages,tau,IC_case,a)
%% ID name generator 
% Generates ID and IDn for an specific simulation using the controling 
% parameters.

%% Read Parameters
% Statistic Used
switch theta
    case {-1}
        statistic = 'BE';
    case {0}
        statistic = 'MB';
    case {1}
        statistic = 'FD';
    otherwise
        error('not a valid theta value')
end
% Elements used
xelements  = num2str(nx);
yelements  = num2str(ny);
vxelements  = num2str(nvx);
vyelements  = num2str(nvy);
% RK stages
RKs = num2str(RK_stages);
% Relaxation time
rtime = num2str(abs(log10(tau)));
% diffusive coef
apha = num2str(a);
% Initial Condition Number
ic = num2str(IC_case);
% Tecplot format
f = '.plt';

%% Generate ID
IDn = [name,'-',statistic,'-',xelements,'x',yelements,'x',vxelements,'x',vxelements,...
    'RK',RKs,'apha',apha,'rtime',rtime,'-','IC',ic,f];
ID = [name,'-',statistic,'-',xelements,'x',yelements,'x',vxelements,'x',vxelements,...
    'RK',RKs,'apha',apha,'rtime',rtime,'-','IC',ic];