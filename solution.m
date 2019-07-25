%function [l_func,i_func,A_out,W] = solution(G,edu,jd,earner,S,eparams,P) 


%% index for parameters

% % estimated
% theta1=eparams(1); % (Kt) efficiency with which investments translate into improvements in cognition
% delta01=eparams(2); % (Kt) natural rate of cognitive decline, by age groups 
% delta02=eparams(3); % (Kt) natural rate of cognitive decline, by age groups 
% delta03=eparams(4); % (Kt) natural rate of cognitive decline, by age groups 
% delta04=eparams(5); % (Kt) natural rate of cognitive decline, by age groups
% theta2=eparams(6); % (Ut) fixed cost of work, in hours - i.e., hours worked per year- estimated in French & Jones
% alpha=eparams(7); % (Ut) the “productive” share of cognitive investment that is not perceived as leisure
% lambda=eparams(8); % (Ut) Weight on leisure (both in TVF and within period UF)
% omega=eparams(9); % (AIME) approximates the ratio of the lowest earnings year to AIME
% nu=eparams(10); % (TVF) coefficient of relative risk aversion, utility
% theta_b=eparms(11); % (TVF) Bequest weight 
% kappa=eparams(12); % (TVF) bequest shifter, in thousands (French & Jones); kappa determines the curvature of the bequest function
% + variance of wage shocks = eparams(13); 
% + variance of health shocks = eparams(14);

% % calibrated
% alpha01=cparams(1); % wage parameters
% alpha02=cparams(2); % wage parameters
% alpha11=cparams(3); % wage parameters
% alpha12=cparams(4); % wage parameters
% alpha2=cparams(5); % wage parameters
% zeta01=cparams(6); % parameters for job demand  
% zeta02=cparams(7); % parameters for job demand  
% zeta03=cparams(8); % parameters for job demand  
% zeta11=cparams(9); % parameters for job demand  
% zeta12=cparams(10); % parameters for job demand  
% tau10=cparams(11); % parameters for job losing probabilities
% tau11=cparams(12); % parameters for job losing probabilities
% tau12=cparams(13); % parameters for job losing probabilities
% tau13=cparams(14); % parameters for job losing probabilities
% tau14=cparams(15); % parameters for job losing probabilities

%% Vector for cognitive investment 
inv_w_vector = linspace(G.inv_min,G.inv_w_max,G.n_coginv);
inv_r_vector = linspace(G.inv_min,G.inv_r_max,G.n_coginv);
            

%% Terminal Value Function
                    
% From bequest function in French 2005 and French & Jones 2011 
TVF = (G.theta_b*(S.SS_A + G.kappa).^(1-G.nu)*G.lambda)/(1-G.nu); % (15*15*5) x 1 = 1125 x 1


%% 1. Loop for time (32 periods, TBD):
tic;
for t = G.n_period-1:-1:1
    t;
    
    % age & age groups
    age = 50 + t;
    agegr = agegr_func(age);

    % Coefficients for Chebyshev Approximation
    if t==G.n_period-1
        Emax = TVF;
            Num= Emax'*S.B; % B= kron(T_A, kron(T_K,T_M)), 1125 x 784
            Den= S.T2; % T2= kron(T2_A, kron(T2_A,T2_A)), 784 x 1
            coeff = Num./Den'; % 1 x 784 
    else
        Emax = W(:,t+1);
            Num= Emax'*S.B; % B= kron(T_A, kron(T_K,T_M))
            Den= S.T2; % T2= kron(T2_A, kron(T2_A,T2_A))
            coeff = Num./Den'; % 1 x 784 
    end

    % 2. Loop over shocks (9 points):
    
    for i = 1:1:G.n_shocks
        i;

        % wage shocks
        shocks_w= S.shocks_w(i);
        shocks_h= S.shocks_h(i);

        % wage (jd-edu-age profile)
        wage=10000; % temporary
         % wage = exp(alpha01 + alpha02*(jd==2) + alpha03*(jd==3) + alpha11*(edu==2) + alpha12*(edu==3) + alpha2*log(1+age) + shocks_w);            

        % 3. Loop over state space (15*15*5 = 1125 points):

        for j = 1:1:(G.n_assets*G.n_cogcap*G.n_aime) % 15*15*5
                j;

            % Individual's assets
              A_j = S.SS_A(j);
              
            % cognitive capital 
              K_j = S.SS_K(j);  
              
            % Individual's aime 
              AIME_j = S.SS_M(j); 
              
            % AIME to PIA 
              PIA_t = PIA_func(AIME_j, age);

            % Labor Income 
              inc_w = wage;
              inc_r = PIA_t;
              
              % vector for cognitive investment 
              inv_min= 0;
              inv_w_max= 24-8-8; % get from data
              inv_r_max= 24-8; % get from data
              inv_w_vector = linspace(inv_min,inv_w_max,G.n_coginv);
              inv_r_vector = linspace(inv_min,inv_r_max,G.n_coginv);
              
              % 4. loop over cognitive investment 

              for k = 1:1:G.n_coginv
                  k; 

                  % cognitive investment by work status - Q. how to combine health investment and job demands for workers?

                   JD = 10; % temporary
                   % JD = zeta01 + zeta02*(jd==2) + zeta03*(jd==3) + zeta11*(edu==2) + zeta12*(edu==3); % TBD 

                   inv_w = inv_w_vector(k) + JD;  
                   inv_r = inv_r_vector(k); 

                   % job probabilities
                     prob_lambda = 0.3; % temporary 
                      % prob_lambda = normcdf(tau10 + tau11*(edu==2) + tau12*(edu==3) + tau13*age + tau14*JD); % probability of losing a job 

                   % Utility by work status: U_t(l_t,K_j)=(T-l_t-a*i_t)^lambda*(K_j)^(1-lambda)
                     u_w(k) = U_func(8, inv_w, K_j, inc_w); %% check 
                     u_r(k) = U_func(0, inv_r, K_j, inc_r); %% check 
                     
                   % transition for cognitive capital: K_(t+1)=?*i_t^total+(1-?_t)K_j+?_t^health
                   % Q. are theta1, delta, shocks same for workers and the retired?
                    K_w_next(k) = K_trans(inv_w, agegr, K_j, shocks_h);
                    K_r_next(k) = K_trans(inv_r, agegr, K_j, shocks_h);
                   
                   % Transition for asset: A_(t+1)=(1+r)A_t+l_t*w_t(JD,Ed,age,?_t^wage)+(1-l_t)PIA(AIME_j,age_t)  
                     A_w_next(k) = A_trans(G, A_j, inc_w, inv_w);
                     A_r_next(k) = A_trans(G, A_j, inc_r, inv_r);
                                                 
                   % Transition for AIME
                     AIME_w_next(k) = AIME_trans(G, age, AIME_j, wage);
                     AIME_r_next(k) = AIME_j; % if retired, AIME doesn't update anymore
                   
                   % polynomial approximation of VF (for all continuous variables)

                     % for workers 
                     Base_w_A=chebpoly_base(S.nA+1, S.d_A*(A_w_next(k) - S.extmin_A) - 1); 
                     Base_w_K=chebpoly_base(S.nK+1, S.d_K*(K_w_next(k) - S.extmin_K) - 1);                                        
                     Base_w_M=chebpoly_base(S.nM+1, S.d_M*(AIME_w_next(k) - S.extmin_M) - 1);
                     Base_w=kron(Base_w_A,kron(Base_w_K,Base_w_M));
                     
                     V_w_next(k) = sum(coeff.*Base_w,2); 
                     
                     % for retired 
                     Base_r_A=chebpoly_base(S.nA+1, S.d_A*(A_r_next(k) - S.extmin_A) - 1); 
                     Base_r_K=chebpoly_base(S.nK+1, S.d_K*(K_r_next(k) - S.extmin_K) - 1);                                        
                     Base_r_M=chebpoly_base(S.nM+1, S.d_M*(AIME_r_next(k) - S.extmin_M) - 1);
                     Base_r=kron(Base_r_A,kron(Base_r_K,Base_r_M));
                     
                     V_r_next(k) = sum(coeff.*Base_r,2); 
                     
                   % validate output (Assets, AIME, Cognitive capital)
                            
                     % for workers           
                     if A_w_next < S.extmin_A
                        V_w_next(k)=NaN;
                     else
                        V_w_next(k)=V_w_next(k);
                     end
                     
                     % for retired          
                     if A_r_next < S.extmin_A
                        V_r_next(k)=NaN;
                     else
                        V_r_next(k)=V_r_next(k);
                     end
                            
                     % Sector-Specific Value Function (for working & retired)
                     V_w(k) = u_w(k) + G.beta * ((1-prob_lambda)*V_w_next(k) + prob_lambda*V_r_next(k));
                     V_r(k) = u_r(k) + G.beta * (V_r_next(k));

                     % Save value function (for work decisions)
                     V_w_aux(k,j,i) = V_w(k);                               
                     V_r_aux(k,j,i) = V_r(k); 
                                    
              end
              
              % optimal investment and max VF
              
              % save optimal 
              [V_w_star, index_w_n] = max(V_w);
              [V_r_star, index_r_n] = max(V_r);
              inv_w_star = inv_w_vector(index_w_n);                            
              inv_r_star = inv_r_vector(index_r_n);
              inv_star_aux = [inv_w_star, inv_r_star]; 
              [V_max, l_index] = max([V_w_star,V_r_star]);
                            
              % save choices
              i_star(j,i,t) = inv_star_aux(l_index);
              l_star(j,i,t) = l_index;
              V_star(j,i,t) = V_max;
              
        end

    end
    
    % Integrate over shocks
    W(:,t) = pi^(-1/2)*V_star(:,:,t)*S.weight;

end
toc