% g^*(x) = r_epsilon(x)
% with Lipschitz constant L_T = max(1/epsilon,1); 
function [T_val, L_T] = T_repsilon(x,epsilon)

         T_val = zeros(size(x));
         x_large = find(abs(x)>epsilon);
         x_small = setdiff(1:length(x), x_large);
         T_val(x_large) = sign(x(x_large));
         T_val(x_small) = x(x_small)/epsilon;
         
         L_T = max(1/epsilon,1);
         
end
