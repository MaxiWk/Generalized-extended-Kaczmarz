fig_folder_name = 'noise_RAc_real';
dir_to_figures = '/Users/maximilianwinkler/Documents/Braunschweig/Forschungsthemen/Stochastic_splitting_methods/Kaczmarz method/Sparse Kaczmarz/ExtendedSparseKaczmarz_octave/tex/figures/';

m = 1000;  % 500 
n = 500;  % 1000
sp = 25;   % 25


num_repeats = 50; 

real_setting = true;

maxiter = 5e5; % Number of iterations % 5e6
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

experiment_description = 'rank-deficient, large noise in R(A) complement';
%noise5_rangeAc

%experiment_description = 'rank deficient, noise split into R(A) and R(A) complement';

%experiment_description = 'rank deficient, only noise in R(A) complement, well conditioned A'; % maxiter = 2*1e5
% with sigma_max = 2, sigma_min = 1 -> very fast with sigma_max = 100,
% convergence, ESREK best, SREK improves very much tue to change in xhat_min, but ESREK surprisingly not

%experiment_description = 'rank deficient, only noise in R(A) complement, well conditioned A and xhat';
% with xhat_max = 6, xhat_min = 1 
% slow convergence, SREK best, ESREK worse than SRK

%experiment_description = 'rank deficient, only noise in R(A) complement, well conditioned A and xhat';  
% with sigma_max = 2, xhat_min = 0, xhat_max = 5 -> very fast
% convergence, ESREK not affected

%experiment_description = 'rank deficient, noise split into R(A) and R(A) complement, well conditioned A and xhat';

%experiment_description = 'dctmatrix, noise split into R(A) and R(A) complement';

%experiment_description = 'rank deficient, consistent, no noise';

%experiment_description = 'rank deficient, consistent, no noise';

disp_instance = false;

stopcrit_sample_pars.length_resAbz_sampled = ceil(m/2);
stopcrit_sample_pars.length_resATz_sampled = ceil(n/2);
stopcrit_sample_pars.min_possible_iter_for_stopping = 4*max(m,n);
T = {T_1,T_2};
L_gstar = [L_gstar_1,L_gstar_2];

data = experiment(n,m,sp,real_setting,lambda_value,T,L_gstar,maxiter,num_repeats,iter_save,rowsamp,colsamp,...
                  writeout,dir_to_figures,fig_folder_name,disp_instance,savestep,stopcrit_sample_pars,method_array,experiment_description);                           

save(fullfile([dir_to_figures '/' fig_folder_name, '/data.mat']), 'data', '-mat');

                              
                              
                              
                              