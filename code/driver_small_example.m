m = 2;
n = 2;
sp = 2;

experiment_description = 'small_example';

num_repeats = 2;

real_setting = false; 

rowsamp = 'rownorms squared';
colsamp = 'colnorms squared';

lambda = 0.5;  

maxiter = 1e4;

T = @(z) z;  % gradient gstar for g(x) = 1/2 ||x||_2^2 + gamma ||x||_1
L_gstar = 1;
T = {T};

writeout = false; 

savestep = 1; 

number_data_points = maxiter;
iter_save = floor(maxiter/number_data_points);

method_array = {'rek'}; 

disp_instance = false;

stopcrit_sample_pars.length_resAbz_sampled = ceil(m/2);
stopcrit_sample_pars.length_resATz_sampled = ceil(n/2);
stopcrit_sample_pars.min_possible_iter_for_stopping = 4*max(m,n);

data = experiment(n,m,sp,real_setting,lambda,T,L_gstar,maxiter,num_repeats,iter_save,rowsamp,colsamp,1,...
                  writeout,disp_instance,savestep,stopcrit_sample_pars,method_array,experiment_description);                           

              