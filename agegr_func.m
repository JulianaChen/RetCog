function [agegr] = agegr_func(age)

    if age <= 58
        agegr=1;
    elseif (59 <= age) && (age <= 66)
        agegr=2;
    elseif (67 <= age) && (age <= 74)
        agegr=3;
    else 
        agegr=4;
    end
    
end