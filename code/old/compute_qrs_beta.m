function beta_upper_bound = compute_qrs_beta(gamma,A)


    sing_values = svd(A);
    sigma_m_sq = min(sing_values(sing_values > eps)); % lowest nonzero sv^2
    
    norm_AF_sq = norm(A,'fro')^2;
    
    r = 4*gamma*sigma_m_sq^2/norm_AF_sq;
    q = 1 + 4*gamma*sigma_m_sq - r;
    s = gamma*norm_AF_sq - 2*gamma*sigma_m_sq + gamma*sigma_m_sq^2/norm_AF_sq + 0.5;
    
    alpha_lower_bound = 4*s - q + sqrt((4*s-q)^2 + 4*r*s - q^2);
    
    beta_upper_bound = 1/sqrt(1+alpha_lower_bound);
    
end


