function [u] = U_func(L_j, inv, K_j, inc)

% current utility is a function of: 
%T: annualized per period endownment of time/leisure- estimated in French & Jones- 4,466 hours 
%l_t: participation in work (0 or 1)- Q. work partipiaction decision itself? or work participation decision * work hours (i.e.theta2 below)?
%theta2: fixed cost of work, in hours - i.e., hours worked per year- estimated in French & Jones- 1,313 hours
%alpha: the “productive” share of cognitive investment that is not perceived as leisure
%inv: time invested in cognitive stimulating activities
%"T- theta2*l_t -?i_t": individual’s time constraint and indicates the quantity of leisure consumed

T= 24-8; % 24-sleep, temporary (4,466 French & Jones 2011) 
%L_j= 8; % 8 hours of work, temporary (1,313 French & Jones 2011)
alpha=0.1; % temporary
lambda1=0.3; % temporary
lambda2=0.3; % temporary
eta=0.5; % fixed savings rate, temporary 
                     
u = lambda1*log(T - L_j- alpha*inv) + lambda2*log(K_j) + (1-lambda1-lambda2)*log((1-eta)*inc);

end

%                      u_w(k) = lambda1*log(T - L_j- alpha*inv_w) + lambda2*log(K_j) + (1-lambda1-lambda2)*log((1-eta)*inc_w); %% check 
%                      u_r(k) = lambda1*log(T - alpha*inv_r) + lambda2*log(K_j) + (1-lambda1-lambda2)*log((1-eta)*inc_r); %% check               
%                      u_w(k)= ((T - L_j- alpha*inv_w).^lambda1) * K_j^(lambda2) + ((1-eta)*inc_w)^(1-lambda1-lambda2); 
%                      u_r(k)= ((T - alpha*inv_r).^lambda1) * K_j^(lambda2) + ((1-eta)*inc_r)^(1-lambda1-lambda2);