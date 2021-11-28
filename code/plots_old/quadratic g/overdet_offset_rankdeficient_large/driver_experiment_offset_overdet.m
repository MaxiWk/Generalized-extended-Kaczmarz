

% Add an offset in the complement of the range of A to the righthand side

% with 'solvable, medium gaussian noise' 
% unfortunately, in both cases the extended method gives more nonzero entries

m = 500;  % 20
n = 200;  % 50
sp = 5;

maxiter = 2*1e5; % Number of iterations
iter_save = floor(maxiter/500);  % each such number of iterations, a data point is added in the error plot

rowsamp = 'rownorms squared';
colsamp = 'colnorms squared';

lambda_value = 10;  

T = @(z) z;  % gradient gstar for g(x) = 1/2 ||x||_2^2 
L = 1;

num_repeats = 60; % 60 Number of repeats over random instances

writeout = false;

savestep = 1; 

method_array = {'rek','srk','srek','esrek'}; 

experiment_description =  'rank-deficient, noise in complement of range';
%experiment_description = 'not solvable, not full rank, noise only in the complement of range';

disp_solution = false;

data = experiment(n,m,sp,lambda_value,T,L,maxiter,num_repeats,iter_save,rowsamp,colsamp,1,...
                              writeout,disp_solution,savestep,method_array,experiment_description);

                              
                              
                              
                              