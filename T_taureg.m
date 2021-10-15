% g^*(x) = r_epsilon(x) + tau/2 ||x||_2^2 
% L_gstar = max(1/epsilon,1) + tau;                         
function T_val = T_taureg(x,epsilon,tau)
  
         T_val = (1./ max(epsilon,abs(x))+tau).* x;   
         
end




