
m = 500;  % 500 
n = 200;  % 200
sp = 5; 


num_repeats = 5; 


real_setting = true;

maxiter = 5*1e6; % Number of iterations
number_data_points = 500;
iter_save = floor(maxiter/number_data_points);  % each such number of iterations, a data point is added in the error plot

rowsamp = 'rownorms squared';
colsamp = 'colnorms squared';

lambda_value = 5;  

T = @(z) z;  % gradient gstar for g(x) = 1/2 ||x||_2^2 + gamma ||x||_1
L_gstar = 1;

writeout = false;

savestep = 1; 

method_array = {'rek','srk','esrek','srek'}; 

%experiment_description = 'rank-deficient, only noise in R(A) complement';
% maxiter = 4*1e6
% Sparsity needs to be high enough: m=60, n=20, sp=5 (4) minimizer is not the chosen vector, but
% is so for sp<4 

%experiment_description = 'rank deficient, noise split into R(A) and R(A) complement';

experiment_description = 'rank deficient, only noise in R(A) complement, well conditioned A'; % maxiter = 2*1e5
% with sigma_max = 2, xhat_min = 1 -> very fast
% convergence, ESREK best, SREK improves very much tue to change in xhat_min, but ESREK surprisingly not

%experiment_description = 'rank deficient, only noise in R(A) complement, well conditioned A and xhat';
% with sigma_max = 100, xhat_min = 1 
% slow convergence, SREK best, ESREK worse than SRK

%experiment_description = 'rank deficient, only noise in R(A) complement, well conditioned A and xhat';  
% with sigma_max = 2, xhat_min = 0, xhat_max = 5 -> very fast
% convergence, ESREK not affected

%experiment_description = 'rank deficient, noise split into R(A) and R(A) complement, well conditioned A and xhat';



disp_instance = false;

stopcrit_sample_pars.length_resAbz_sampled = ceil(m/2);
stopcrit_sample_pars.length_resATz_sampled = ceil(n/2);
stopcrit_sample_pars.min_possible_iter_for_stopping = 4*max(m,n);

data = experiment(n,m,sp,real_setting,lambda_value,T,L_gstar,maxiter,num_repeats,iter_save,rowsamp,colsamp,1,...
                  writeout,disp_instance,savestep,stopcrit_sample_pars,method_array,experiment_description);
                              

                              
                              
                              
                              