fig_folder_name ='impulsive_noise_complex_m1000n500_unif_lambda5_medium_conditioned_larger_noise_uniform_p';
dir_to_figures = '/Users/maximilianwinkler/Documents/Braunschweig/Forschungsthemen/Stochastic_splitting_methods/Kaczmarz method/Sparse Kaczmarz/ExtendedSparseKaczmarz_octave/tex/figures/';

m = 1000; % DCT matrix: m=2000  % otherwise: 500/200
n = 500; % DCT matrix: n=1000  % otherwise: 200/500
sp = 25; 
%sp = min(5,n); % as before

real_setting = true; 

maxiter = 2e7;  

num_repeats = 10; 

writeout = true;

iter_save = ceil(maxiter/500);  % each such number of iterations, a data point is added in the error plot

rowsamp = 'rownorms squared';
colsamp = 'colnorms squared'; 



% f(x) = epsilon/2 ||x||_2^2 + gamma ||x||_1  Nonsparse: Set gamma=0.
% T = grad g^*, see below for different T
lambda = 5;   
%gamma = 0.1;
epsilon = 1e-4;    % 0.01
tau = 0.001;    
%mu = 1;
r_epsilon_factor = 0.1;
shrinkage_par = 10; 

T_1 = @(x) x;
L_gstar_1 = 1;

T_2 = @(x) T_taureg(x,epsilon,tau);  
L_gstar_2 = 1/epsilon + tau;

T = {T_1, T_2};
L_gstar = [L_gstar_1, L_gstar_2];

%T = @(x) x; % L_gstar = 1;
%L_gstar = 1;
%T = @(x) max(abs(x)-gamma,0).*sign(x); % L_gstar = 1;
%L_gstar = 1;

%T = @(x) T_repsilon(x,epsilon); % L_gstar = max(1/epsilon,1);
%L_gstar = max(1/epsilon,1);
%T = @(x) T_mureg(x,epsilon,mu);  % L_gstar = mu* max(1/epsilon,1) + 1;
%L_gstar = mu* max(1/epsilon,1) + 1;
%T = @(x) T_repsilon_shrinkage(x,epsilon,r_epsilon_factor,shrinkage_par);  % L_gstar = r_epsilon_factor* max(1/epsilon,1) + 1;

savestep = 1; 

method_array = {'rek','srk','grek_1','grek_2'}; 

%experiment_description = 'impulsive noise, rank deficient, uniform, well conditioned';
experiment_description = 'impulsive noise, rank deficient, uniform, medium conditioned';
%experiment_description = 'impulsive noise, rank deficient, uniform, bad conditioned';

median_res = zeros(maxiter,length(method_array));






disp_instance = false;

stopcrit_sample_pars.length_resAbz_sampled = ceil(m/2);
stopcrit_sample_pars.length_resATz_sampled = ceil(n/2);
stopcrit_sample_pars.min_possible_iter_for_stopping = 4*max(m,n);

data = experiment(n,m,sp,real_setting,lambda,T,L_gstar,maxiter,num_repeats,iter_save,rowsamp,colsamp,...
                          writeout,dir_to_figures,fig_folder_name,disp_instance,savestep,stopcrit_sample_pars,method_array,experiment_description);

save(fullfile([dir_to_figures '/' fig_folder_name '/data.mat']), 'data', '-mat');
                                             
                   
                           



















