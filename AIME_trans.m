function [AIME_next] = AIME_trans(G, age, AIME_j, wage)

omega_t=0.1; % temporary; approximates the ratio of the lowest earnings year to AIME
             % estimate by simulating wage (not earnings) histories with the model developed in French 2005
                            
    if age <= 55
        AIME_next= (1+G.g)*AIME_j + 1/35*wage; % omega_t=0 for workers aged 55 and younger
    elseif (55 < age) && (age <= 60)
        AIME_next= (1+G.g)*AIME_j + 1/35*max(0, wage - omega_t*(1+G.g)*AIME_j);
    else
        AIME_next= (1+G.g)*AIME_j + 1/35*max(0, wage - omega_t*AIME_j); % earnings accrued after age 60 are not rescaled
    end

end