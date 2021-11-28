% g^*(x) = mu r_epsilon(x) + ||x||_2^2
% L_gstar = mu* max(1/epsilon,1) + 1;
function T_val = T_mureg(x,epsilon,mu)
         T_val = zeros(size(x));
         x_large = find(abs(x)>epsilon);
         x_small = setdiff(1:length(x), x_large);
         T_val(x_large) = mu*sign(x(x_large)) + x(x_large);
         T_val(x_small) = (mu+epsilon)/epsilon* x(x_small);         
end