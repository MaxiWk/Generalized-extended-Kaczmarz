% g^*(x) = r_epsilon_factor* r_epsilon(x) + ||S_lambda(x)||_2^2
% L_gstar = r_epsilon_factor* max(1/epsilon,1) + 1;
function T_val = T_repsilon_shrinkage(x,epsilon,r_epsilon_factor,shrinkage_par)
  
         T_val = zeros(size(x));
         x_large = find(abs(x)>epsilon);
         x_small = x_small = setdiff(1:length(x), x_large);
         T_val(x_large) = sign(x(x_large));
         T_val(x_small) = x(x_small)/epsilon;
         
         T_val = r_epsilon_factor* T_val + max(abs(x - shrinkage_par),0).* sign(x);
         
end
