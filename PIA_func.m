function [PIA_t] = PIA_func(AIME_j, age)

%PIA
PB1 = 612*12; % bend point1, annualized, as of 2004
PB2 = 1689*12; % bend point2, annualized, as of 2004 

    if AIME_j < PB1
        PIA_t = 0.9*AIME_j;
    elseif (PB1 <= AIME_j) && (AIME_j < PB2)
        PIA_t = 0.9*PB1 + 0.32*(AIME_j-PB1);
    else
        PIA_t = 0.9*PB1 + 0.32*(PB2-PB1) + 0.15*(AIME_j-PB2);
    end

% Age adjusted PIA:
% If age < 62 -> PIA = 0 (can't claim)
% If age = 62 -> 75%
% If age = 63 -> 80%
% If age = 64 -> 86.7%
% If age = 65 -> 93.3%
% If age > 65 -> 100%

    if (age == 65) && (age < 62)
        PIA_t = 0.933*PIA_t;
    elseif age == 64
        PIA_t = 0.867*PIA_t;
    elseif age == 63
        PIA_t = 0.8*PIA_t;
    elseif age == 62
        PIA_t = 0.75*PIA_t;
    else
        PIA_t = 0;
    end

end