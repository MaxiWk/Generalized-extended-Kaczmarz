

% Add an offset in the complement of the range of A to the righthand side

% with 'solvable, medium gaussian noise' 
% unfortunately, in both cases the extended method gives more nonzero entries


m = 5;  
n = 10;
sp = 2;

maxiter = 1e4; % Number of iterations
num_rand_repeats = 1; 
iter_save = floor(maxiter/500);  % each such number of iterations, a data point is added in the error plot

rowsamp = 'rownorms squared';
colsamp = 'colnorms squared'; 

% lambda: f(x) = epsilon/2 ||x||_2^2 + gamma ||x||_1  Nonsparse: Set gamma=0.
lambda = 5;    % lambda = 5

% gamma: g(x) = epsilon/2 ||x||_2^2 + gamma ||x||_1  Nonsparse g: Set gamma=0.
gamma = 0.1;    
T = @(z) max(abs(z)-gamma,0).*sign(z);  % gradient gstar for g(x) = 1/2 ||x||_2^2 + gamma ||x||_1
% T = @(z) z;
L_gstar = 1;


num_repeats = 5; % 60 Number of repeats over random instances

writeout = false; 

savestep = 1; 

method_array = {'rek','srk','srek'}; 

experiment_description = 'rank-deficient, gaussian noise in complement of range and then uniform noise added'; %'rank deficient, small gaussian noise';

median_res = zeros(maxiter,length(method_array));

xhat_list_srek = zeros(n, num_rand_repeats, num_repeats);
yhat_list_srek = zeros(m, num_rand_repeats, num_repeats);
xhat_list_srk = zeros(n, num_rand_repeats, num_repeats);

for rand_repeats = 1:num_rand_repeats
  
  seeds_indices = 1:num_rand_repeats;
  
  if rand_repeats == 1
    disp_instance = 1;
  else 
    disp_instance = 0;
  end

  data = experiment(n,m,sp,lambda,T,L_gstar,maxiter,num_repeats,iter_save,rowsamp,colsamp,seeds_indices(rand_repeats),...
                              writeout,disp_instance,savestep,method_array,experiment_description);
                        
  xhat_list = data.xhat_list;
  yhat_list = data.yhat_list;  
  xhat_list_srek(:,rand_repeats,:) = xhat_list(:,:,1); 
  yhat_list_srek(:,rand_repeats,:) = yhat_list(:,:,1); 
  xhat_list_srk(:,rand_repeats,:) =  xhat_list(:,:,2);   
  
end
                   
                   
                              