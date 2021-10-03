% requires for all i: (xplus)_i \neq 0 
function [lambda_list, x_list] = first_steps_homotopy_method(A,b,num_timesteps)
  
  Aplus = pinv(A);
  xplus = Aplus*b;
  s = sign(xplus);
  P = eye(size(A,2)) - Aplus*A;
  
  Ps = P*sign(xplus);
  sign_change = find( Ps~= 0);
  t_max = min( abs(xplus(sign_change))./abs(Ps(sign_change)) );
  
  lambda_list = linspace(0,t_max,num_timesteps);
  
  x_list = zeros(length(xplus),num_timesteps);
  x_list(:,1) = xplus;
  
  for i = 0:num_timesteps-1      
      x_list(:,i+1) = xplus - i/(num_timesteps-1)*t_max*Ps;
  end
  
end

