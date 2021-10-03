

% Add an offset in the complement of the range of A to the righthand side

% with 'solvable, medium gaussian noise' 
% unfortunately, in both cases the extended method gives more nonzero entries


m = 200;  
n = 500;
sp = min(5,n);

maxiter = 1e6; % Number of iterations
num_rand_repeats = 1; 
%iter_save = 10;
iter_save = ceil(maxiter/500);  % each such number of iterations, a data point is added in the error plot

rowsamp = 'rownorms squared';
colsamp = 'colnorms squared'; 



% f(x) = epsilon/2 ||x||_2^2 + gamma ||x||_1  Nonsparse: Set gamma=0.
% T = grad g^*, see below for different T
lambda = 10;    % lambda = 5
gamma = 0.1;
epsilon = 0.1;    % >0
tau = 0.1;    % gamma=0.05 besser als gamma=0 für m=500,n=200 - echt?!
mu = 1;
r_epsilon_factor = 0.1;
shrinkage_par = 1;
%T = @(x) T_taureg(x,epsilon,tau);  % L_gstar = max(1/epsilon,1) + tau; 
%L_gstar = max(1/epsilon,1) + tau; 
%T = @(x) x; % L_gstar = 1;
%L_gstar = 1;
%T = @(x) max(abs(x)-gamma,0).*sign(x); % L_gstar = 1;




%T = @(x) T_repsilon(x,epsilon); % L_gstar = max(1/epsilon,1);
%L_gstar = max(1/epsilon,1);
T = @(x) T_mureg(x,epsilon,mu);  % L_gstar = mu* max(1/epsilon,1) + 1;
L_gstar = mu* max(1/epsilon,1) + 1; 
%T = @(x) T_repsilon_shrinkage(x,epsilon,r_epsilon_factor,shrinkage_par);  % L_gstar = r_epsilon_factor* max(1/epsilon,1) + 1;

num_repeats = 5; % 60 Number of repeats over random instances

writeout = false;

savestep = 1; 

%method_array = {'rek','srk','srek'}; 
method_array = {'srk','rek','srek'}; 

experiment_description = 'impulsive noise'; %'rank deficient, small gaussian noise';

median_res = zeros(maxiter,length(method_array));

xhat_srek = zeros(n,num_repeats,num_rand_repeats);
yhat_srek = zeros(m,num_repeats,num_rand_repeats);






for rand_repeats = 1:num_rand_repeats
  
  seeds_indices = 1:num_rand_repeats;

  data = experiment(n,m,sp,lambda,T,L_gstar,maxiter,num_repeats,iter_save,rowsamp,colsamp,seeds_indices(rand_repeats),...
                              writeout,savestep,method_array,experiment_description);

  xhat_list = data.xhat_list;                            
  xhat_srek(:,:,rand_repeats) = xhat_list(:,:,2);    
  yhat_list= data.yhat_list;
  yhat_srek(:,:,rand_repeats) = yhat_list(:,:,2);   
  
end
                   
                   
                           



















