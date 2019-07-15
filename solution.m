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


%% Terminal Value Function

% age at TVF
age_TVF = 82; % initial age in 2004 (51-70) -> final age in 2016 (63-82), Final age: TBD

% assets
assets = S.SS_A; 
                    
% TVF (bequest function, French 2005, French & Jones 2011) 
TVF = (G.theta_b*(assets + G.kappa).^(1-G.nu)*G.lambda)/(1-G.nu); % 1125 (15*15*5) x 1


%% 1. Loop for time (32 periods, TBD):

for t = G.n_period-1:-1:1
    
    % age
    age = 50 + t;
    
    % age groups 
    if age <= 58
        agegr=1;
    elseif 59 <= age <= 66
        agegr=2;
    elseif 67 <= age <= 74
        agegr=3;
    else 
        agegr=4;
    end 
     
    t;

    % Coefficients for Chebyshev Approximation
    if t==G.n_period-1
        Emax = TVF;
            Num= Emax'*S.B; % B= kron(T_A, kron(T_K,T_M))
            Den= T2; % T2= kron(T2_A, kron(T2_A,T2_A))
            coeff = Num./Den'; % 1 x 784 
    else
        Emax = W(:,:,t+1);
            Num= Emax'*S.B; % B= kron(T_A, kron(T_K,T_M))
            Den= T2; % T2= kron(T2_A, kron(T2_A,T2_A))
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
                

        % 3. Loop over assets (15 points):

        for j = 1:1:(G.n_assets*G.n_aime*G.n_cogcap) % 15*15*5
                j;

            % Individual's assets
              A_j = S.SS_A(j);
              
            % Individual's aime 
              AIME_j = S.SS_M(j); 
              
            % AIME to PIA 
              PB1= 612*12; % bend point1, annualized, as of 2004
              PB2= 1689*12; % bend point2, annualized, as of 2004 

              if AIME_j < PB1 
                 PIA_t = 0.9*AIME_j;
              elseif PB1 <= AIME_j < PB2
                 PIA_t = 0.9*PB1 + 0.32*(AIME_j-PB1);
              else 
                 PIA_t = 0.9*PB1 + 0.32*(PB2-PB1) + 0.15*(AIME_j-PB2);
              end
              
            % cognitive capital 
              K_j = S.SS_K(j);  

              % vector for cognitive investment 
              inv_min= 0;
              inv_max= 168 ; % get from data
              inv_vector = linspace(inv_min,inv_max,G.n_coginv);
                      
              % 4. loop over cognitive investment 

              for k = 1:1:G.n_coginv
                  k; 

                  % cognitive investment by work status - Q. how to combine health investment and job demands for workers?

                   JD = 10; % temporary
                   % JD = zeta01 + zeta02*(jd==2) + zeta03*(jd==3) + zeta11*(edu==2) + zeta12*(edu==3); % TBD 

                   inv_w = inv_vector(k) + JD;  
                   inv_r = inv_vector(k); 

                   % transition for cognitive capital: K_(t+1)=?*i_t^total+(1-?_t)K_j+?_t^health
                   % Q. are theta1, delta, shocks same for workers and the retired?

                     theta1=1; % efficiency with which investments translate into improvements in cognition, temporary
                     delta01=0.4; % age-related natural depreciation rate for age 51-58, temporary
                     delta02=0.3; % age-related natural depreciation rate for age 59-66, temporary
                     delta03=0.2; % age-related natural depreciation rate for age 67-74, temporary
                     delta04=0.1; % age-related natural depreciation rate for age 75-82, temporary

                     K_w_next = theta1*inv_w + ...
                              + (1-delta01)*(agegr==1)*K_j...
                              + (1-delta02)*(agegr==2)*K_j...
                              + (1-delta03)*(agegr==3)*K_j...
                              + (1-delta04)*(agegr==4)*K_j...
                              + shocks_h; % estimate delta by age groups (4 deltas); 
                          
                     K_r_next = theta1*inv_r + ...
                              + (1-delta01)*(agegr==1)*K_j...
                              + (1-delta02)*(agegr==2)*K_j...
                              + (1-delta03)*(agegr==3)*K_j...
                              + (1-delta04)*(agegr==4)*K_j...
                              + shocks_h; % estimate delta by age groups (4 deltas); 

                   % job probabilities

                     prob_lamba = 0.3; % temporary 
                      % prob_lamba = normcdf(tau10 + tau11*(edu==2) + tau12*(edu==3) + tau13*age + tau14*JD); % probability of losing a job 

                   % Utility by work status: U_t(l_t,K_j)=(T-l_t-a*i_t)^lambda*(K_j)^(1-lambda)
                            
                     % current utility is a function of: 
                     %T: annualized per period endownment of time/leisure- estimated in French & Jones- 4,466 hours 
                     %l_t: participation in work (0 or 1)- Q. work partipiaction decision itself? or work participation decision * work hours (i.e.theta2 below)?
                     %theta2: fixed cost of work, in hours - i.e., hours worked per year- estimated in French & Jones- 1,313 hours
                     %alpha: the “productive” share of cognitive investment that is not perceived as leisure
                     %inv: time invested in cognitive stimulating activities
                     %"T- theta2*l_t -?i_t": individual’s time constraint and indicates the quantity of leisure consumed

                     T= 4466; % temporary (French & Jones 2011) 
                     L_j= 1313; % temporary (French & Jones 2011)
                     alpha=0.1; % temporary
                     lambda=0.5; % temporary
                     eta=0.5; % fixed savings rate, temporary 

                     inc_w = wage;
                     inc_r = PIA_t;
                     
                     u_w(k) = lambda*log(T - L_j- alpha*inv_w) + (1-lambda)*log(K_j); % + (1-eta)*inc_w; %% check 
                     u_r(k) = lambda*log(T - alpha*inv_r) + (1-lambda)*log(K_j); % + (1-eta)*inc_r; %% check 
                     
%                      u_w(k)= ((T - L_j- alpha*inv_w).^lambda) * K_j^(1-lambda) + (1-eta)*inc_w; 
%                      u_r(k)= ((T - alpha*inv_r).^lambda) * K_j^(1-lambda) + (1-eta)*inc_r;
                   
                   % Transition for asset: A_(t+1)=(1+r)A_t+l_t*w_t(JD,Ed,age,?_t^wage)+(1-l_t)PIA(AIME_j,age_t)
                     
                     iota=10; % price of health investment, temporary 
                     
                     A_w_next = (1+G.r)*A_j + eta*inc_w - iota*inv_w; % add fixed ratio of savings, minus health inv cost?
                     A_r_next = (1+G.r)*A_j + eta*inc_r - iota*inv_r; % add fixed ratio of savings, minus health inv cost?
                            
                   % Transition for AIME
                    
                     omega_t=0.1; % temporary; approximates the ratio of the lowest earnings year to AIME
                                  % estimate by simulating wage (not earnings) histories with the model developed in French 2005
                            
                     if age <= 55 
                     AIME_w_next= (1+G.g)*AIME_j + 1/35*wage; % omega_t=0 for workers aged 55 and younger
                     elseif 55 < age <= 60
                     AIME_w_next= (1+G.g)*AIME_j + 1/35*max(0, wage - omega_t*(1+G.g)*AIME_j);
                     else 
                     AIME_w_next= (1+G.g)*AIME_j + 1/35*max(0, wage - omega_t*AIME_j); % earnings accrued after age 60 are not rescaled 
                     end
                   
                     AIME_r_next = AIME_j; % if retired, AIME doesn't update anymore
                   
                   % polynomial approximation of VF (for all continuous variables)

                     % for workers 
                     Base_w_A=chebpoly_base(S.nA+1, S.d_A*(A_w_next - S.extmin_A) - 1); 
                     Base_w_K=chebpoly_base(S.nK+1, S.d_K*(K_w_next - S.extmin_K) - 1);                                        
                     Base_w_M=chebpoly_base(S.nM+1, S.d_M*(AIME_w_next - S.extmin_M) - 1);
                     Base_w=kron(Base_w_A,kron(Base_w_K,Base_w_M));
                     
                     V_w_next = sum(coeff.*Base_w,2); 
                     
                     % for retired 
                     Base_r_A=chebpoly_base(S.nA+1, S.d_A*(A_r_next - S.extmin_A) - 1); 
                     Base_r_K=chebpoly_base(S.nK+1, S.d_K*(K_r_next - S.extmin_K) - 1);                                        
                     Base_r_M=chebpoly_base(S.nM+1, S.d_M*(AIME_r_next - S.extmin_M) - 1);
                     Base_r=kron(Base_r_A,kron(Base_r_K,Base_r_M));
                     
                     V_r_next = sum(coeff.*Base_r,2); 
                     
                   % save output (Assets, AIME, Cognitive capital )
                            
                     % for workers 
                     A_w_next(k)=A_w_next;
                               
                     if A_w_next < S.extmin_A
                        V_w_next(k)=NaN;
                     else
                        V_w_next(k)=V_w_next;
                     end
                     
                     AIME_w_next(k)=AIME_w_next;
                     K_w_next(k)=K_w_next;
                     
                     % for retired 
                     A_r_next(k)=A_r_next;
                               
                     if A_r_next < S.extmin_A
                        V_r_next(k)=NaN;
                     else
                        V_r_next(k)=V_r_next;
                     end

                     AIME_r_next(k)=AIME_r_next;
                     K_r_next(k)=K_r_next;
                            
                     % Sector-Specific Value Function (for working & retired)
                     V_w(k) = u_w(k) + G.beta * ((1-prob_lamba)*V_w_next(k) + prob_lamba*V_r_next(k));
                     V_r(k) = u_r(k) + G.beta * (V_r_next(k));

                     % Save value function (for work decisions)
                     V_w_aux(k,j,i) = V_w(k);                               
                     V_r_aux(k,j,i) = V_r(k); 
                                    
              end
        end
                         
                     % optimal investment and max VF
                       % save optimal 
                       [V_w_star, index_w_n] = max(V_w);
                       [V_r_star, index_r_n] = max(V_r);
                       inv_w_star = inv_vector(index_w_n);                            
                       inv_r_star = inv_vector(index_r_n);
                       inv_star_aux = [inv_w_star, inv_r_star]; 
                       [V_star, l_index] = max([V_w_star,V_r_star]);
                            
                       % save choices
                       i_star(k,j,i,t) = inv_star_aux(l_index);
                       l_star(k,j,i,t) = l_index;
                       V_star(k,j,i,t) = V_star;

                    end

                end

        
    % Integrate over shocks
    W(:,l,t) = pi^(-1/2)*V_star(:,:,l,t)*S.weight;
  
    % save policy functions 
                        