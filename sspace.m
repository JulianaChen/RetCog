% function [S] = sspace(eparams,G)

%% Shocks - Q. Persistent shocks? how to differntiate wage shocks & health shocks?

sigma_w = 0.43; % wage shocks
sigma_h = 0.43; % health shocks

[e, wt] = GaussHermite(G.Ne);
eps_w = sqrt(2)*e*sigma_w; % error vector
eps_h = sqrt(2)*e*sigma_h; % error vector

% vcv = [sigma_eps^2,0;0,sigma_v^2]; % DON'T NEED THIS NOW, THEY'RE INDEPENDENT
% detV = det(vcv);
% detV = det(vcv);
% detR = det(R);

shocks_w = kron(eps_w,ones(length(eps_w),1)); % 9x1
shocks_h = kron(eps_h,ones(length(eps_h),1)); % 9x1

weight = kron(wt,wt); %kron(wt, kron(wt,wt)); % 9x1

% Basis for Income Shocks
zeps_w= 2*(eps_w-eps_w(1))/(eps_w(G.Ne,1)-eps_w(1))-1; 
zeps_h= 2*(eps_h-eps_h(1))/(eps_h(G.Ne,1)-eps_h(1))-1; 

Teps_w=chebpoly_base(G.Ne-1,zeps_w);
Teps_h=chebpoly_base(G.Ne-1,zeps_h);

T2eps_w = diag(Teps_w'*Teps_w);
T2eps_h = diag(Teps_h'*Teps_h);

%% Discrete Variables

work = [1 0];
SS_L = work; 

%% Continuous Variables (CHEVYSHEV GRID)

% assets_lb  = 0; 
% assets_int = 211000; % 50%
% assets_ub  = 545000; % 75% 
% assets1 = linspace(assets_lb,assets_int,G.n_assets-3);
% assets2 = linspace(assets_int,assets_ub,4);
% assets = [assets1,assets2([2:4])];
%
% cogcap_lb  = 0; 
% cogcap_int = 17; % 50%
% cogcap_ub  = 27; % 100%
% cogcap = linspace(cogcap_lb,cogcap_ub,G.n_cogcap);
% 
% aime_lb  = 0; 
% aime_int = 5000; % 50%, check SSA data
% aime_ub  = 20000; % 100%, check SSA data
% aime = linspace(aime_lb,aime_ub,G.n_aime);

% Testing:
asset1=0;
asset2=545000;
cogcap1=0;
cogcap2=27; 
aime1=0;
aime2=20000;

% Chevyshev Approximation
[assets,nA,extmin_A,extmax_A,d_A,T_A,T2_A] = cheby_values(G.n_assets,asset2,asset1);
[cogcap,nK,extmin_K,extmax_K,d_K,T_K,T2_K] = cheby_values(G.n_cogcap,cogcap2,cogcap1);
[aime,nM,extmin_M,extmax_M,d_M,T_M,T2_M] = cheby_values(15,aime2,aime1); %% aime has different size - should all continuous variables be the same size?

%% SS for chevyshev - Q. how does this work? 

SS_A = assets;
SS_K = cogcap;
SS_M = aime;

% Expand Vector (to calculate function)
% SS_A = repmat(assets',[length(cogcap)*length(aime) 1]);
% SS_K = repmat(kron(cogcap',ones(length(aime),1)),[length(assets) 1]);
% SS_M = kron(aime',ones([length(assets)*length(cogcap),1]));

%% output

S = struct(...
    'eps_w',eps_w,'eps_h',eps_h,'Teps_w',Teps_w,'Teps_h',Teps_h,'T2eps_w',T2eps_w,'T2eps_h',T2eps_h,...
    'shocks_w',shocks_w,'shocks_h',shocks_h,'weight',weight,...
    'extmin_A',extmin_A,'extmax_A',extmax_A,'extmin_K',extmin_K,'extmax_K',extmax_K,...
    'extmin_M',extmin_M,'extmax_M',extmax_M,'nA',nA,'nK',nK,'nM',nM,...
    'd_A', d_A, 'd_K', d_K, 'd_M', d_M, 'assets', assets, 'cogcap', cogcap, 'aime', aime,...
    'SS_A', SS_A, 'SS_K', SS_K, 'SS_M', SS_M,'T_A',T_A,'T_K',T_K,'T_M',T_M,...
    'T2_A',T2_A,'T2_K',T2_K,'T2_M',T2_M,'SS_L',SS_L);

% end 