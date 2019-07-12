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
% alpha02=cparams(2);
% alpha11=cparams(3);
% alpha12=cparams(4);
% alpha2=cparams(5);
% zeta01=cparams(6); % parameters for job demand  
% zeta02=cparams(7);
% zeta03=cparams(8);
% zeta11=cparams(9);
% zeta12=cparams(10);
% tau10=cparams(11); % parameters for job losing probabilities
% tau11=cparams(12);
% tau12=cparams(13);
% tau13=cparams(14);
% tau14=cparams(15);


%% Initial AIME & PIA 

% Initial AIME (draw from joint distribution of relevant variables in SSA earnings history) 

AIME_t= 10000; % approximately average of 35 years of earnings 

% Initial PIA (monthly pension benefit amount)
                                    
PB1= 612*12; % annualized, as of 2004
PB2= 1689*12; % annualized, as of 2004 
                                        
if AIME_t < PB1 
PIA_t = 0.9*AIME_t;
elseif PB1 <= AIME_t < PB2
PIA_t = 0.9*PB1 + 0.32*(AIME_t-PB1);
else 
PIA_t = 0.9*PB1 + 0.32*(PB2-PB1) + 0.15*(AIME_t-PB2);
end

%% Terminal Value Function

% age at TVF

age_TVF = 82; % initial age in 2004 (51-70) -> final age in 2016 (63-82), Final age: TBD

% assets

assets = S.SS_A; 
                    
% TVF (bequest function, French 2005, French & Jones 2011) - Q. do we use this formula for bequest function?

% TVF = repmat((G.theta_b*(assets + G.kappa).^(1-G.nu)*G.lambda)'/(1-G.nu),1,2); %15 (assets) x 2 (work)? 
TVF = repmat((G.theta_b*(assets + G.kappa).^(1-G.nu)*G.lambda)/(1-G.nu),1,2); % 3375 x 2 


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
    
%     % Coefficients for Chebyshev Approximation
%     if t==G.n_period-1
%         Emax = TVF;
%         for x = 1:1:(G.n_work) % 
%             Num_A(x,:) = Emax(:,x)'*S.T_A; % Ts (polynomial base for two dimensional approximation) are all same
%             Den_A = S.T2_A;
%             coeff_A(x,:) = Num_A(x,:)./Den_A';
%             Num_K(x,:) = Emax(:,x)'*S.T_K;
%             Den_K = S.T2_K;
%             coeff_K(x,:) = Num_K(x,:)./Den_K';
%             Num_M(x,:) = Emax(:,x)'*S.T_M; 
%             Den_M = S.T2_M;
%             coeff_M(x,:) = Num_M(x,:)./Den_M';
%         end
%     else
%         Emax = W(:,:,t+1);
%         for x = 1:1:(G.n_work)
%             Num_A(x,:) = Emax(:,x)'*S.T_A;
%             Den_A = S.T2_A;
%             coeff(x,:) = Num_A(x,:)./Den_A';
%             Num_K(x,:) = Emax(:,x)'*S.T_K;
%             Den_K = S.T2_K;
%             coeff(x,:) = Num_K(x,:)./Den_K';
%             Num_M(x,:) = Emax(:,x)'*S.T_M; 
%             Den_M = S.T2_M;
%             coeff(x,:) = Num_M(x,:)./Den_M';
%         end
%     end


    % Coefficients for Chebyshev Approximation
    if t==G.n_period-1
        Emax = TVF;
        for x = 1:1:(G.n_work) 
            Num_A(x,:) = Emax(:,x)'*S.B_A; % T_A -> B_A
            Den_A = kron(T2_A, kron(T2_A,T2_A)); % S.T2_A -> kron(T2_A, kron(T2_A,T2_A))
            coeff_A(x,:) = Num_A(x,:)./Den_A';
            Num_K(x,:) = Emax(:,x)'*S.B_K;
            Den_K = kron(T2_K, kron(T2_K,T2_K));
            coeff_K(x,:) = Num_K(x,:)./Den_K';
            Num_M(x,:) = Emax(:,x)'*S.B_M; 
            Den_M = kron(T2_M, kron(T2_M,T2_M));
            coeff_M(x,:) = Num_M(x,:)./Den_M';
        end
    else
        Emax = W(:,:,t+1);
        for x = 1:1:(G.n_work)
           Num_A(x,:) = Emax(:,x)'*S.B_A; 
            Den_A = kron(T2_A, kron(T2_A,T2_A));
            coeff_A(x,:) = Num_A(x,:)./Den_A';
            Num_K(x,:) = Emax(:,x)'*S.B_K;
            Den_K = kron(T2_K, kron(T2_K,T2_K));
            coeff_K(x,:) = Num_K(x,:)./Den_K';
            Num_M(x,:) = Emax(:,x)'*S.B_M; 
            Den_M = kron(T2_M, kron(T2_M,T2_M));
            coeff_M(x,:) = Num_M(x,:)./Den_M';
        end
    end

    coeff = coeff_A; % temporary 

    % 2. Loop over assets (15 points):

        for i = 1:1:G.n_assets
                i;

        % Individual's assets
          A_j = S.SS_A(i);
        
  
            % 3. loop for shocks (9 points):
                for j = 1:1:G.n_shocks
                        j;
                    
                % wage shocks
                    shocks_w= S.shocks_w(j);
                    shocks_h= S.shocks_h(j);

                % wage (jd-edu-age profile)
                
%                   wage = exp(alpha01 + alpha02*(jd==2) + alpha03*(jd==3) + alpha11*(edu==2) + alpha12*(edu==3) + alpha2*log(1+age) + shocks_w);
                    wage=10000; % temporary 
                    
                    
                      % 4. loop for work status (2 points):
                        for x = 1:1:G.n_work
                                x; 
                                
                        % current state discrete variables (work):
                          L_j = S.SS_L(x);
                         
                          
                            % 5. loop over cognitive capital (15 points)
                              for k = 1:1:G.n_cogcap
                                      k;
                                      
                              % cognitive capital 
                                K_t = S.SS_K(k);  
                                
                              % vector for cognitive investment 
                                inv_min= 0;
                                inv_max= 168 ; % get from data
                                inv_vector = linspace(inv_min,inv_max,G.n_coginv);
                                
                                
                                % 6. loop over cognitive investment 
                                
                                  for m = 1:1:G.n_coginv
                                      m; 
                                      
                                    % cognitive investment by work status - Q. how to combine health investment and job demands for workers?

%                                       JD = zeta01 + zeta02*(jd==2) + zeta03*(jd==3) + zeta11*(edu==2) + zeta12*(edu==3); % TBD 
                                        JD = 10; % temporary
             
                                        inv_w = inv_vector(m) + JD;  
                                        inv_r = inv_vector(m); 
                                      
                                    % transition for cognitive capital: K_(t+1)=?*i_t^total+(1-?_t)K_t+?_t^health
                                    
                                      % Q. are theta1, delta, shocks same for workers and the retired?
                                            
                                        theta1=1; % efficiency with which investments translate into improvements in cognition, temporary
                                        delta01=0.4; % age-related natural depreciation rate for age 51-58, temporary
                                        delta02=0.3; % age-related natural depreciation rate for age 59-66, temporary
                                        delta03=0.2; % age-related natural depreciation rate for age 67-74, temporary
                                        delta04=0.1; % age-related natural depreciation rate for age 75-82, temporary

                                        K_next = (theta1*inv_w)*L_j + (theta1*inv_r)*(1-L_j) + ...
                                                + (1-delta01)*(agegr==1)*K_t...
                                                + (1-delta02)*(agegr==2)*K_t...
                                                + (1-delta03)*(agegr==3)*K_t...
                                                + (1-delta04)*(agegr==4)*K_t...
                                                + shocks_h; % estimate delta by age groups (4 deltas); 
                                            
                                    % job probabilities
                                        
%                                       prob_lamba = normcdf(tau10 + tau11*(edu==2) + tau12*(edu==3) + tau13*age + tau14*JD); % probability of losing a job
                                        prob_lamba = 0.3; % temporary 
                                        
                                    % Utility by work status: U_t(l_t,K_t)=(T-l_t-?i_t )^?*(K_t)^(1-?)
                                      % current utility is a function of: 
                                       %T: annualized per period endownment of time/leisure- estimated in French & Jones- 4,466 hours Q.do we estimate T? 
                                       %l_t: participation in work (0 or 1)- Q. work partipiaction decision itself? or work participation decision * work hours (i.e.theta2 below)?
                                       %theta2: fixed cost of work, in hours - i.e., hours worked per year- estimated in French & Jones- 1,313 hours
                                       %alpha: the “productive” share of cognitive investment that is not perceived as leisure
                                       %inv: time invested in cognitive stimulating activities
                                       %"T- theta2*l_t -?i_t": individual’s time constraint and indicates the quantity of leisure consumed
                                        
                                        T= 4466; % temporary (French & Jones 2011) 
                                        theta2= 1313; % temporary (French & Jones 2011)
                                        alpha=0.1; % temporary
                                        lambda=0.5; % temporary
                                        eta=0.5; % fixed savings rate, temporary 
                                        iota=10; % price of health investment, temporary 

                                        inc = L_j*wage + (1-L_j)*PIA_t;
                                        
                                        u_w(m)= ((T - theta2*L_j- alpha*inv_w).^lambda) * K_t^(1-lambda) + (1-eta)*inc - iota*inv_w; 
                                        u_r(m)= ((T - theta2*(1-L_j)- alpha*inv_r).^lambda) * K_t^(1-lambda) + (1-eta)*inc - iota*inv_r;
                                        
                                    % Transition for asset: A_(t+1)=(1+r)A_t+l_t*w_t(JD,Ed,age,?_t^wage)+(1-l_t)PIA(AIME_t,age_t)
                                    
                                        A_next = (1+G.r)*A_j + eta*inc; % add fixed ratio of savings, deduct health inv cost?
                                    
                                    % For workers
                                        
                                        % Transition for AIME
                                        
                                        omega_t=0.1; % temporary; approximates the ratio of the lowest earnings year to AIME
                                                     % estimate by simulating wage (not earnings) histories with the model developed in French 2005
                                        
                                        if age <= 55
                                        AIME_next= (1+G.g)*AIME_t + 1/35*wage; % omega_t=0 for workers aged 55 and younger
                                        elseif 55 < age <= 60
                                        AIME_next= (1+G.g)*AIME_t + 1/35*max(0, wage - omega_t*(1+G.g)*AIME_t);
                                        else 
                                        AIME_next= (1+G.g)*AIME_t + 1/35*max(0, wage - omega_t*AIME_t); % earnings accrued after age 60 are not rescaled 
                                        end
                                        
                                        % polynomial approximation of VF (for all continuous variables)
%                                     
                                        Base_A=chebpoly_base(S.nA+1, S.d_A*(A_next - S.extmin_A) - 1);
                                        Base_K=chebpoly_base(S.nK+1, S.d_K*(K_next - S.extmin_K) - 1);                                        
                                        Base_M=chebpoly_base(S.nK+1, S.d_K*(K_next - S.extmin_K) - 1);
                                      
                                        Base=kron(Base_A, kron(Base_K, Base_M));
                                        V_w_next=sum(coeff.*Base,2);
                                        
                                        % polynomial approximation of VF (for all continuous variables)- using cheby_approx
                                        V_w_next=cheby_approx(coeff,13,S.extmin_A,S.extmin_K,S.extmin_M,S.d_A,S.d_K,S.d_M,A_next,K_next,AIME_next);
                                        
                                        % save output (Assets, AIME, Cognitive capital )

                                        A_w_next(k)=A_next;
                                        if A_next < S.extmin_A
                                            V_w_next(k)=NaN;
                                        else
                                            V_w_next(k)=V_w_next;
                                        end
                                        
                                        AIME_w_next(k)=AIME_next;
                                        
                                        K_w_next(k)=K_next;
                                      

                                    % For retired

                                        % Transition for AIME: if retired, AIME doesn't update anymore

                                        % polynomial approximation of VF (for all continuous variables)
%                                         Base_A=chebpoly_base(S.nA+1, S.d_A*(A_next - S.extmin_A) - 1);
%                                         Base_K=chebpoly_base(S.nK+1, S.d_K*(K_next - S.extmin_K) - 1);                                        
%                                         Base_M=chebpoly_base(S.nM+1, S.d_M*(M_next - S.extmin_M) - 1);
                                        
%                                         V_r_next_A = sum(coeff(x_next,:).*Base_A,2); %cheby_approx
%                                         V_r_next_K = sum(coeff(x_next,:).*Base_K,2); %cheby_approx
%                                         V_r_next_M = sum(coeff(x_next,:).*Base_M,2); %cheby_approx

                                        % save output (AIME & Assets)

                                        AIME_r_next(k)=AIME_next;

                                        A_r_next(k)=A_next;
                                        if A_next < S.extmin_A
                                            V_r_next(k)=NaN;
                                        else
                                            V_r_next(k)=V_r_next;
                                        end
                                        
                                        K_r_next(k)=K_next;
                                        
                                    % Sector-Specific Value Function (for working & retired)
                                        V_w(k) = u_w(k) + G.beta * (prob_lamba*V_w_next(k) + (1-prob_lamba)*V_r_next(k));
                                        V_r(k) = u_r(k) + G.beta * (V_r_next(k));
                                        
%                                     % Save value function (for work decisions)
%                                         V_w_aux(l,k,x,i,j) = V_w(k);                               
%                                         V_r_aux(l,k,x,i,j) = V_r(k); 
                                            
                                  end 
                                                
                              end 
                        
                              % Choose cognitive investment 
                             
                        end 
                    
                        % Choose work decisions 
                
                    % optimal investment and max VF
            
                    % save choices
            
            end 
            
            % integrate over shocks 
            
    end 
    
    % save policy functions 
    
end     
