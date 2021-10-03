

% Add an offset in the complement of the range of A to the righthand side

% with 'solvable, medium gaussian noise' 
% unfortunately, in both cases the extended method gives more nonzero entries

m = 500;
n = 200;
sp = 5;

maxiter = 3*1e5; % Number of iterations
iter_save = floor(maxiter/500);  % each such number of iterations, a data point is added in the error plot

rowsamp = 'rownorms squared';
colsamp = 'colnorms squared'; 

lambda_value = 5;  

gamma = 0.1;    % 0.1
T = @(z) max(abs(z)-gamma,0).*sign(z);  % gradient gstar for g(x) = 1/2 ||x||_2^2 + gamma ||x||_1
L_gstar = 1;

num_repeats = 3; % 60 Number of repeats over random instances

writeout = false;

savestep = 1; 

method_array = {'rek','srk','srek'}; 

experiment_description = 'medium uniform noise';

data = experiment(n,m,sp,lambda_value,T,L_gstar,maxiter,num_repeats,iter_save,rowsamp,colsamp,1,...
                              writeout,savestep,method_array,experiment_description);

                       
                             
                   
                   
                              