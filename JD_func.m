function [JD] = JD_func(cparams, jd, edu)

zeta01=cparams(6); % parameters for job demand  
zeta02=cparams(7); % parameters for job demand  
zeta03=cparams(8); % parameters for job demand  
zeta11=cparams(9); % parameters for job demand  
zeta12=cparams(10); % parameters for job demand

JD = zeta01 + zeta02*(jd==2) + zeta03*(jd==3) + zeta11*(edu==2) + zeta12*(edu==3); % TBD

end