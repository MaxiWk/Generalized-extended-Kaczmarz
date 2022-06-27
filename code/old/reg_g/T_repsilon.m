% g^*(x) = r_epsilon(x)
% L_gstar = max(1/epsilon,1); 
function T_val = T_repsilon(x,epsilon)
         T_val = zeros(size(x));
         x_large = find(abs(x)>epsilon);
         x_small = x_small = setdiff(1:length(x), x_large);
         T_val(x_large) = sign(x(x_large));
         T_val(x_small) = x(x_small)/epsilon;
end
