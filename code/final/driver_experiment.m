% driver for experiment (i) 
% Maximilian Winkler, maximilian.winkler@tu-bs.de

dir_to_figures = '';
fig_folder_name = '';

m = 10;  % 500 
n = 20;  % 1000
sp = 5;   % 25

num_repeats = 2; 

real = false;

maxiter = 2e5; % iterations
number_data_points = 500;
iter_save = floor(maxiter/number_data_points);  % each such number of iterations, a data point is added in the error plot

rowsamp = 'uniform';
colsamp = 'uniform';

lambda_value = 5;  

T_1 = @(z) z;  % gradient gstar for g(x) = 1/2 ||x||_2^2 
L_gstar_1 = 1;

% T for Frank's example
epsilon = 0.01;
tau = 0.001;
T_2 = @(x) T_taureg(x,epsilon,tau);
L_gstar_2 = 1/epsilon + tau;

writeout = true; 

savestep = 1; 

method_array = {'rek', 'srk', 'grek_1'}; 

%experiment_description = 'rank-deficient, medium noise in R(A) complement';

experiment_description = 'rank-deficient, large noise in R(A) complement';

stopcrit_sample_pars.length_resAbz_sampled = ceil(m/2);
stopcrit_sample_pars.length_resATz_sampled = ceil(n/2);
stopcrit_sample_pars.min_possible_iter_for_stopping = 4*max(m,n);
T = {T_1,T_2};
L_gstar = [L_gstar_1,L_gstar_2];

data = experiment(n, m, sp, real, lambda_value, T, L_gstar, maxiter, num_repeats, iter_save, rowsamp, colsamp,...
                  writeout, dir_to_figures, fig_folder_name, savestep, stopcrit_sample_pars, ...
                  method_array,experiment_description);                           

save('data.mat', 'data', '-mat');

                              
                              
                              
                              