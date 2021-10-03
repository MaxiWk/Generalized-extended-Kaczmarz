

% Add an offset in the complement of the range of A to the righthand side

% with 'solvable, medium gaussian noise' 
% unfortunately, in both cases the extended method gives more nonzero entries


m = 500;  
n = 200;
sp = 5;

maxiter = 5*1e5; % Number of iterations
num_rand_repeats = 1; 
iter_save = ceil(maxiter/500);  % each such number of iterations, a data point is added in the error plot

rowsamp = 'rownorms squared';
colsamp = 'colnorms squared'; 

% lambda: f(x) = epsilon/2 ||x||_2^2 + gamma ||x||_1  Nonsparse: Set gamma=0.
lambda = 5;    % lambda = 5

gamma = 0.1;    
T = @(z) max(abs(z)-gamma,0).*sign(z);  % gradient gstar for g(x) = 1/2 ||x||_2^2 + gamma ||x||_1
L_gstar = 1;

num_repeats = 5; % 60 Number of repeats over random instances

writeout = false;

savestep = 1; 

method_array = {'rek','srk','srek'}; 

experiment_description = 'not solvable, noise in complement of range and gaussian noise added';

median_res = zeros(maxiter,length(method_array));

xhat_srek = zeros(n,num_repeats,num_rand_repeats);
yhat_srek = zeros(m,num_repeats,num_rand_repeats);

for rand_repeats = 1:num_rand_repeats
  
  seeds_indices = 1:num_rand_repeats;

  data = experiment(n,m,sp,lambda,T,L_gstar,maxiter,num_repeats,iter_save,rowsamp,colsamp,seeds_indices(rand_repeats),...
                              writeout,savestep,method_array,experiment_description);
%{
  xhat_list = data.xhat_list;                            
  xhat_srek(:,:,rand_repeats) = xhat_list(:,:,2);    
  yhat_list= data.yhat_list;
  yhat_srek(:,:,rand_repeats) = yhat_list(:,:,2);   
%}
end
                   
                   
                              