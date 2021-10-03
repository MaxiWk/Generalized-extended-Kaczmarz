

% Add an offset in the complement of the range of A to the righthand side

% with 'solvable, medium gaussian noise' 
% unfortunately, in both cases the extended method gives more nonzero entries

m = 800;  % 20
n = 600;  % 50
sp = ceil(n/20); % 5

num_repeats = 60; 

real_setting = true;

maxiter = 2*1e5; % Number of iterations
iter_save = floor(maxiter/500);  % each such number of iterations, a data point is added in the error plot

rowsamp = 'rownorms squared';
colsamp = 'colnorms squared';

lambda_value = 10;  

T = @(z) z;  % gradient gstar for g(x) = 1/2 ||x||_2^2 + gamma ||x||_1
L_gstar = 1;

writeout = false;

savestep = 1; 

method_array = {'rek','srk','esrek','srek'}; 

%experiment_description =  'rank-deficient, noise in complement of range';

experiment_description = 'not solvable, noise only in the complement of range, well conditioned A';
% with sigma_min = 1, sigma_max = 2 -> very fast convergence, ESREK best

%experiment_description = 'not solvable, noise in complement of the range,
%well conditioned A and xhat';  % with sigma_max = 2, xhat_min = 1 -> very fast
%convergence, ESREK best, SREK improves very much tue to change in xhat_min, but ESREK surprisingly not

%experiment_description = 'not solvable, noise in complement of the range,
%well conditioned A and xhat';  % with sigma_max = 100, xhat_min = 1 
% slow convergence, SREK best, ESREK worse than SRK

%experiment_description = 'not solvable, noise in the complement of range,
%well conditioned A and xhat';  % with sigma_max = 2, xhat_min = 0, xhat_max = 5 -> very fast
%convergence, ESREK not affected

%experiment_description = 'not solvable, noise in the complement of range and noise in range, well conditioned A and xhat';

disp_instance = false;

data = experiment(n,m,sp,real_setting,lambda_value,T,L_gstar,maxiter,num_repeats,iter_save,rowsamp,colsamp,1,...
                              writeout,disp_instance,savestep,method_array,experiment_description);

                              
                              
                              
                              