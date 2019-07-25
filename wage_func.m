function [wage] = wage_func(cparams, jd, edu, age, shocks_w)

alpha01=cparams(1); % wage parameters
alpha02=cparams(2); % wage parameters
alpha11=cparams(3); % wage parameters
alpha12=cparams(4); % wage parameters
alpha2=cparams(5); % wage parameters

wage = exp(alpha01 + alpha02*(jd==2) + alpha03*(jd==3) + alpha11*(edu==2) + alpha12*(edu==3) + alpha2*log(1+age) + shocks_w);            

end