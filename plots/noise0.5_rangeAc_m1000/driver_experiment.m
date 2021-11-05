dir_to_folder_with_figures = 'plots/noise0.5_rangeAc_m1000';

m = 1000;  % 500 
n = 500;  % 200
sp = 50;   % 20


num_repeats = 50;  % 5

real_setting = true;

maxiter = 1e7; % Number of iterations
number_data_points = 500;
iter_save = floor(maxiter/number_data_points);  % each such number of iterations, a data point is added in the error plot

rowsamp = 'rownorms squared';
colsamp = 'colnorms squared';

lambda_value = 5;  

T_1 = @(z) z;  % gradient gstar for g(x) = 1/2 ||x||_2^2 + gamma ||x||_1
L_gstar_1 = 1;

T = {T_1};
L_gstar = {L_gstar_1};

writeout = false; 

savestep = 1; 

method_array = {'rek','srk','grek_1'}; 

experiment_description = 'rank-deficient, medium noise in R(A) complement';
% noise0.5_rangeAc

%experiment_description = 'rank-deficient, large noise in R(A) complement';
% noise5_rangeAc

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


disp_instance = false;

stopcrit_sample_pars.length_resAbz_sampled = ceil(m/2);
stopcrit_sample_pars.length_resATz_sampled = ceil(n/2);
stopcrit_sample_pars.min_possible_iter_for_stopping = 4*max(m,n);

data = experiment(n,m,sp,real_setting,lambda_value,T,L_gstar,maxiter,num_repeats,iter_save,rowsamp,colsamp,1,...
                  writeout,disp_instance,savestep,stopcrit_sample_pars,method_array,experiment_description);                           

save(fullfile(dir_to_folder_with_figures, 'data.mat'), 'data', '-mat');
save_figures(dir_to_folder_with_figures)
                              
                              
                              
                              