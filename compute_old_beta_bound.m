function old_beta_bound = compute_old_beta_bound(A,lambda,xhat)

    gamma = compute_gamma(A,lambda,xhat);
    
    frac = 1/(1+2*gamma*norm(A,'fro')^2);
    
    old_beta_bound = sqrt( 0.5* (1 - sqrt(1-frac)) );

end