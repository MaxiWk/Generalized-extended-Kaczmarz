

soft = @(x) max(abs(x)-lambda,0).*sign(x); % Soft thresholding 
grad_Huber = @(x) (1./ max(epsilon,abs(x))+tau).* x;   

% GERK-(a,d) method (sparse in x, L2-proj of b onto R(A))
grad_fstar = soft;
L_fstar = 1;
grad_gstar = @(x) x;
L_gstar = 1;
gerk_ad = GERK_method('gerk_ad', grad_fstar, L_fstar, grad_gstar, L_gstar);

% GERK-(b,d) method (sparse in x, approximated L1-proj of b onto R(A))
grad_fstar = soft;
L_fstar = 1;
grad_gstar = grad_Huber;
L_gstar = 1/epsilon + tau;

gerk_bd = GERK_method('gerk_bd', grad_fstar, L_fstar, grad_gstar, L_gstar);
