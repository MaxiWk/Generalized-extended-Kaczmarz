
% Add an offset in the complement of the range of A to the righthand side

% with 'solvable, medium gaussian noise' 
% unfortunately, in both cases the extended method gives more nonzero entries

m = 2000;  % 20
n = 1000;  % 50
sp = ceil(n/20); % 5

num_repeats = 2; 

real_setting = true;

maxiter = 1e5; % Number of iterations
number_data_points = 500;
iter_save = floor(maxiter/number_data_points);  % each such number of iterations, a data point is added in the error plot

rowsamp = 'rownorms squared';
colsamp = 'colnorms squared';

lambda_value = 10;  

T = @(z) z;  % gradient gstar for g(x) = 1/2 ||x||_2^2 + gamma ||x||_1
L_gstar = 1;

writeout = false;

savestep = 1; 

%method_array = {'rek','srk','esrek','srek'};
method_array = {'esrek','srek'};

experiment_description = 'eyedct, noise in R(A)c and a bit of noise in R(A)';

disp_instance = false;

stopcrit_sample_pars.length_resA_sampled = ceil(m/2);
stopcrit_sample_pars.length_resAT_sampled = ceil(n/2);
stopcrit_sample_pars.min_possible_iter_for_stopping = 4*max(m,n);

data = experiment(n,m,sp,real_setting,lambda_value,T,L_gstar,maxiter,num_repeats,iter_save,rowsamp,colsamp,1,...
                  writeout,disp_instance,savestep,stopcrit_sample_pars,method_array,experiment_description);
                              

                              
                              
                              
                              