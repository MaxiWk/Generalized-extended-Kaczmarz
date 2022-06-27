

% Add an offset in the complement of the range of A to the righthand side

% with 'solvable, medium gaussian noise' 
% unfortunately, in both cases the extended method gives more nonzero entries


m = 500;  % Verhalten dreht sich um, wenn m und n getauscht werden -??
n = 200;
sp = 5;

maxiter = 5*1e5; % Number of iterations
num_rand_repeats = 1; 
iter_save = floor(maxiter/500);  % each such number of iterations, a data point is added in the error plot

rowsamp = 'rownorms squared';
colsamp = 'colnorms squared'; 

% lambda: f(x) = epsilon/2 ||x||_2^2 + gamma ||x||_1  Nonsparse: Set gamma=0.
lambda = 5;    % lambda = 5
% gamma: g(x) = epsilon/2 ||x||_2^2 + gamma ||x||_1  Nonsparse g: Set gamma=0.
gamma = 0.05;    % gamma=0.05 besser als gamma=0 für m=500,n=200 - echt?!

num_repeats = 5; % 60 Number of repeats over random instances

writeout = false;

savestep = 1; 

%method_array = {'rek','srk','srek'}; 
method_array = {'srk','rek','srek'}; 

experiment_description = 'rank-deficient, gaussian noise in complement of range and then uniform noise added'; %'rank deficient, small gaussian noise';

median_res = zeros(maxiter,length(method_array));

xhat_srek = zeros(n,num_repeats,num_rand_repeats);
yhat_srek = zeros(m,num_repeats,num_rand_repeats);

for rand_repeats = 1:num_rand_repeats
  
  seeds_indices = 1:num_rand_repeats;

  data = experiment(n,m,sp,lambda,gamma,maxiter,num_repeats,iter_save,rowsamp,colsamp,seeds_indices(rand_repeats),...
                              writeout,savestep,method_array,experiment_description);

  xhat_list = data.xhat_list;                            
  xhat_srek(:,:,rand_repeats) = xhat_list(:,:,2);    
  yhat_list= data.yhat_list;
  yhat_srek(:,:,rand_repeats) = yhat_list(:,:,2);   
  
end
                   
                   
                              