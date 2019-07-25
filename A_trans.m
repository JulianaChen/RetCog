function [A_next] = A_trans(G, A_j, inc, inv)

eta=0.5; % fixed savings rate, temporary 
iota=10; % price of investments, temporary

A_next = (1+G.r)*A_j + eta*inc - iota*inv;

end