% g^*(x) = r_epsilon(x) + tau/2 ||x||_2^2 
% with Lipschitz constant L_T = 1/epsilon + tau;                         
function [T_val, L_T] = T_taureg(x,epsilon,tau)
  
         T_val = (1./ max(epsilon,abs(x))+tau).* x;   
         
         L_T = 1/epsilon + tau;
         
end




