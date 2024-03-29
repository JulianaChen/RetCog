function [K_next] = K_trans(inv, agegr, K_j, shocks_h)

theta1=1; % efficiency with which investments translate into improvements in cognition, temporary
delta01=0.4; % age-related natural depreciation rate for age 51-58, temporary
delta02=0.3; % age-related natural depreciation rate for age 59-66, temporary
delta03=0.2; % age-related natural depreciation rate for age 67-74, temporary
delta04=0.1; % age-related natural depreciation rate for age 75-82, temporary

K_next = theta1*inv + ...
    + (1-delta01)*(agegr==1)*K_j...
    + (1-delta02)*(agegr==2)*K_j...
    + (1-delta03)*(agegr==3)*K_j...
    + (1-delta04)*(agegr==4)*K_j...
    + shocks_h;

end
