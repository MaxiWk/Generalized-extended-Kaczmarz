

% Add an offset in the complement of the range of A to the righthand side

% with 'solvable, medium gaussian noise' 
% unfortunately, in both cases the extended method gives more nonzero entries

m = 3;  
n = 2;  
sp = 2;

maxiter = 1e5; % Number of iterations 
time_per_repeat = 5; % time to run if 'experiment_time'
iter_save = ceil(maxiter/100);  % each such number of iterations, a data point is added in the error plot
iter_save = 100;

rowsamp = 'rownorms squared';
colsamp = 'colnorms squared'; 

lambda_value = 0.1;  
gamma = 0.1;

%T = @(z) max(abs(z)-gamma,0).*sign(z);  % gradient gstar for g(x) = 1/2 ||x||_2^2 + gamma ||x||_1
T = @(z) z;
L_gstar = 1;

num_repeats = 1; % 60 Number of repeats over random instances

writeout = false;
disp_instance = false;

savestep = 1; 

method_array = {'srk'} %{'lin_breg','srk','srek'}; %{'srek','det_sek'}; 

experiment_description = 'rank-deficient, noise in complement of range'; %'rank deficient, medium uniform noise';

experiment_description = 'not solvable, medium gaussian noise';

data = experiment(n,m,sp,lambda_value,T,L_gstar,maxiter,num_repeats,iter_save,rowsamp,colsamp,1,...
                              writeout,disp_instance,savestep,method_array,experiment_description);



% run each algorithm for x seconds                          
%data = experiment_time(n,m,sp,lambda_value,T,L_gstar,maxiter,time_per_repeat,num_repeats,iter_save,rowsamp,colsamp,1,...
%                              writeout,disp_instance,savestep,method_array,experiment_description);
                              
                              
                              
                              
                              
                              
                              
                              
                              
         

         