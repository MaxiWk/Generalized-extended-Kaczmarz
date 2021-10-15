
m = 2000; % DCT matrix: m=2000  % otherwise: 500/200
n = 1000; % DCT matrix: n=1000  % otherwise: 200/500
sp = 50; % Frank's examples: ceil(m/20) 
%sp = min(5,n); % as before


real_setting = true;

maxiter = 4*1e6;  % DCT matrix: 4*1e6, otherwise 2*1e5


num_repeats = 2; % 60 Number of repeats over random instances 

%iter_save = 10;
iter_save = ceil(maxiter/500);  % each such number of iterations, a data point is added in the error plot

rowsamp = 'rownorms squared';
colsamp = 'colnorms squared'; 



% f(x) = epsilon/2 ||x||_2^2 + gamma ||x||_1  Nonsparse: Set gamma=0.
% T = grad g^*, see below for different T
lambda = 5;    % lambda = 5
%gamma = 0.1;
epsilon = 0.01;    % >0 
tau = 0.001;    
%mu = 1;
r_epsilon_factor = 0.1;
shrinkage_par = 1; 


T = @(x) T_taureg(x,epsilon,tau);  
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

method_array = {'srek','esrek'}; 

experiment_description = 'dctmatrix, noise split into R(A) and R(A) complement'; %'impulsive noise 10%, full rank'; 

median_res = zeros(maxiter,length(method_array));






disp_instance = false;

stopcrit_sample_pars.length_resAbz_sampled = ceil(m/2);
stopcrit_sample_pars.length_resATz_sampled = ceil(n/2);
stopcrit_sample_pars.min_possible_iter_for_stopping = 4*max(m,n);

data = experiment(n,m,sp,real_setting,lambda,T,L_gstar,maxiter,num_repeats,iter_save,rowsamp,colsamp,1,...
                          writeout,disp_instance,savestep,stopcrit_sample_pars,method_array,experiment_description);

               
                   
                           



















