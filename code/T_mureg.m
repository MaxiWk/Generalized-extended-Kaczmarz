% g^*(x) = mu r_epsilon(x) + ||x||_2^2
% with Lipschitz constant L_T = mu* max(1/epsilon,1) + 1;
function [T_val, L_T] = T_mureg(x,epsilon,mu)

         T_val = zeros(size(x));
         x_large = find(abs(x)>epsilon);
         x_small = ~ismembc(T_val,x_large);
         T_val(x_large) = mu*sign(x(x_large)) + x(x_large);
         T_val(x_small) = (mu+epsilon)/epsilon* x(x_small);      
         
         L_T = mu* max(1/epsilon,1) + 1;
         
end