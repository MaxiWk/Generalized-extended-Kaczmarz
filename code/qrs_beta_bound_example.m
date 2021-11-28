rng(0)

for i = 1:5
    A = randn(10,20);
    sp = 5;

    xhat = sparserandn(size(A,2),sp);
    lambda = 1;

    gamma = compute_gamma(A,lambda,xhat)
    beta_upper_bound = compute_qrs_beta(gamma,A)
    old_beta_bound = compute_old_beta_bound(A,lambda,xhat)
end