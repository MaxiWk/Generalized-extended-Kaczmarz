function gamma = compute_gamma(A,lambda,xhat)

    abs_xhat = abs(xhat);
    xhat_min = min(abs_xhat(abs_xhat>eps));
    xhat_min_factor = (xhat_min + 2*lambda)/ xhat_min;
    
    % compute sigma expression
    sigma_tilde_sqr = inf;
    n = size(A,2);
    for num_cols = 1:n
        combos = nchoosek(1:n,num_cols);
        for combo_nr = 1:size(combos,1)
            A_small = A(:,combos(combo_nr,:));
            sing_values = svd(A_small);
            min_sigma_sqr = min(sing_values(sing_values>eps))^2;
            sigma_tilde_sqr = min(sigma_tilde_sqr, min_sigma_sqr);
        end
    end
       
    gamma = 1/sigma_tilde_sqr* xhat_min_factor;

end