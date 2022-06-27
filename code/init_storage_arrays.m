% Initialize storage arrays for results
num_methods = length(method_ids); 
[res, lsres, grad_zfunctional, err_to_sparse, err_to_moorepi, nonzero_entries]...
                     = deal( zeros(floor(maxiter/iter_save), num_repeats, num_methods) );
err_to_sparse_stopped = zeros(num_repeats, num_methods);   
err_to_sparse_end_div_stopped = zeros(num_repeats, num_methods); 
x_last_test_instance = zeros(n,num_methods);
[length_resAbz_sampled, length_resATz_sampled] = deal(stopcrit_sample_pars.length_resAbz_sampled,...
                                                             stopcrit_sample_pars.length_resATz_sampled);
resAbz_sampled = zeros(1, stopcrit_sample_pars.length_resAbz_sampled);
resATz_sampled = zeros(1, stopcrit_sample_pars.length_resATz_sampled);
[resAbz_mean, resATz_mean] = deal( zeros(1,maxiter) );
iterstop_list = zeros(num_repeats,num_methods);
xstop_list = zeros(n, num_repeats, num_methods);
[resAbz_list, resATz_list] = deal( zeros(floor(maxiter/iter_save), num_repeats, num_methods) );
[resAbz_mean_list, resATz_mean_list] = deal( zeros(floor(maxiter/iter_save), num_repeats, num_methods) );
time = zeros(num_repeats, num_methods);

min_rel_cond_A = inf;
max_rel_cond_A = -inf;
