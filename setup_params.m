%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% MATLAB CODE RETIREMENT AND COGNITION %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Gauss-Hermite Points
Ne=3;

%Discount rate
beta=0.7;

%Interest rate
r=0.07;

%Growth rate of average earnings in the overall economy 
g=0.026; % for the period of 2004-2016 

%Parameters for TVF (from French & Jones 2011, French 2005) - calibrate 
nu = 7.49; % coefficient of relative risk aversion, utility (French & Jones)
kappa = 0.444; % bequest shifter, in thousands (French & Jones); kappa determines the curvature of the bequest function
lambda = 0.649; % consumption weight 
theta_b = 0.0223; % bequest weight

%State parameters
n_shocks = 9; % income shocks, health shocks 
n_period = 32; % age 51-70 in 2004 -> age 63 - 82 in 2016?
n_pop = 3000; % TBD
n_assets = 15; % 15 points of support
n_aime = 5; % 5 points of support 
n_coginv= 10; % 10 points of support
n_cogcap= 15; % 15 points of support
n_work= 2; % 2 points of support

% Health investment 
inv_min= 0;
inv_w_max= 24-8-8; % get from data
inv_r_max= 24-8; % get from data

%simulation parameters
Eps=randn(3,n_pop,n_period);

G = struct('Ne',Ne,'beta',beta,'r',r,'g',g,'Eps',Eps,...
    'nu',nu,'kappa',kappa,'lambda',lambda,'theta_b',theta_b,...
    'n_period',n_period,'n_shocks',n_shocks,'n_pop',n_pop,'n_work',n_work,...
    'n_assets',n_assets,'n_aime',n_aime,'n_coginv',n_coginv,'n_cogcap',n_cogcap,...
    'inv_min',inv_min,'inv_w_max',inv_w_max,'inv_r_max',inv_r_max);

%% Estimated Parameters (import from excel)

% [eparams, enames] = xlsread(paramfile, 'eParams');

%% Calibrated Parameters (import from excel)

% [cparams, cnames] = xlsread(paramfile, 'cParams');

% for k = 1:length(cparams)
%      CurrVarname=cell2mat(cnames(k));
%      CurrValue=cparams(k);
%      P.(CurrVarname)=CurrValue;
% end

%% Initial Conditions

edu_levels = [1:3];
jd_levels = [1:3];
earner_type = [1 2];

edu_levels_ext= [repmat(1,3);repmat(2,3);repmat(3,3)];
jd_levels_ext = repmat([1:3]',3);

types = [kron(earner_type',ones(length(edu_levels_ext),1)) ...
        [edu_levels_ext([1:9],1);edu_levels_ext([1:9],1)] ... 
        repmat(edu_levels',[length(jd_levels)*length(earner_type) 1])];
    
    