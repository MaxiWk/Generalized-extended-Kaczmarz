
m = 100; %20;  
n = 50; %50;
sp = min(5,n); 

maxiter = 5*1e6; %5*1e4; % Number of iterations

num_repeats = 20; % 60 Number of repeats over random instances 

num_rand_repeats = 1; % how many times with different seeds 
%iter_save = 10;
iter_save = ceil(maxiter/500);  % each such number of iterations, a data point is added in the error plot

rowsamp = 'rownorms squared';
colsamp = 'colnorms squared'; 



% f(x) = epsilon/2 ||x||_2^2 + gamma ||x||_1  Nonsparse: Set gamma=0.
% T = grad g^*, see below for different T
lambda = 5;    % lambda = 5
%gamma = 0.1;
epsilon = 0.0001;    % >0 
tau = 1;    
%mu = 1;
r_epsilon_factor = 0.1;
shrinkage_par = 1; 


T = @(x) T_taureg(x,epsilon,tau);  % L_gstar = max(1/epsilon,1) + tau; 
L_gstar = max(1/epsilon,1) + tau; 
%T = @(x) x; % L_gstar = 1;
%L_gstar = 1;
%T = @(x) max(abs(x)-gamma,0).*sign(x); % L_gstar = 1;
%L_gstar = 1;

%T = @(x) T_repsilon(x,epsilon); % L_gstar = max(1/epsilon,1);
%L_gstar = max(1/epsilon,1);
%T = @(x) T_mureg(x,epsilon,mu);  % L_gstar = mu* max(1/epsilon,1) + 1;
%L_gstar = mu* max(1/epsilon,1) + 1;
%T = @(x) T_repsilon_shrinkage(x,epsilon,r_epsilon_factor,shrinkage_par);  % L_gstar = r_epsilon_factor* max(1/epsilon,1) + 1;

writeout = false;

savestep = 1; 

method_array = {'rek','srk','srek'}; 

experiment_description = 'large impulsive noise'; %'rank deficient, small gaussian noise';

median_res = zeros(maxiter,length(method_array));

xhat_srek = zeros(n,num_repeats,num_rand_repeats);
yhat_srek = zeros(m,num_repeats,num_rand_repeats);






for rand_repeats = 1:num_rand_repeats
  
  seeds_indices = 1:num_rand_repeats;
  
  disp_instance = 0;

  data = experiment(n,m,sp,lambda,T,L_gstar,maxiter,num_repeats,iter_save,rowsamp,colsamp,seeds_indices(rand_repeats),...
                              writeout,disp_instance,savestep,method_array,experiment_description);

  xhat_list = data.xhat_list;                            
  xhat_srek(:,:,rand_repeats) = xhat_list(:,:,2);    
  yhat_list= data.yhat_list;
  yhat_srek(:,:,rand_repeats) = yhat_list(:,:,2);   
  
end
                   
                   
                           



















