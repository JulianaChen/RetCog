function [prob_lambda] = prob_func(cparams, edu, age, JD)

tau10=cparams(11); % parameters for job losing probabilities
tau11=cparams(12); % parameters for job losing probabilities
tau12=cparams(13); % parameters for job losing probabilities
tau13=cparams(14); % parameters for job losing probabilities
tau14=cparams(15); % parameters for job losing probabilities

prob_lambda = normcdf(tau10 + tau11*(edu==2) + tau12*(edu==3) + tau13*age + tau14*JD); % probability of losing a job

end