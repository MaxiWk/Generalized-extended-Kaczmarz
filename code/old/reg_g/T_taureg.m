% g^*(x) = r_epsilon(x) + tau/2 ||x||_2^2 
% L_gstar = max(1/epsilon,1) + tau;                         
function T_val = T_taureg(x,epsilon,tau)
         T_val = zeros(size(x));
         
         x_large = find(abs(x)>1+epsilon*tau);
         T_val(ismembc(T_val,x_large)) = sign(x(x_large)) + tau*x(x_large);
         T_val(~ismembc(T_val,x_large)) = (1/epsilon + tau)*x(x_small);
         
         %x_large = find(abs(x)>1+epsilon*tau);
         %x_small = find(abs(x)<=1+epsilon*tau);
         %T_val(x_large) = sign(x(x_large)) + tau*x(x_large);
         %T_val(x_small) = (1/epsilon + tau)*x(x_small);
end

