% driver for experiment (i) 
% Maximilian Winkler, maximilian.winkler@tu-bs.de
% Lionel Ngoupeyou Tondji, l.ngoupeyou-tondji.tu-bs.de

method_ids = {'rek','srk','gerk_ad','gerk_bd'}; 

m = 1000;   % 1000
n = 500;  % 500
sp = min(m,n)/20; 

real = false; 

num_repeats = 50; % 50
maxiter = 2e5; %2e5 

number_data_points = 500;
iter_save = floor(maxiter/number_data_points);  % each such number of iterations, a data point is added in the error plot

rowsamp = 'uniform';
colsamp = 'uniform'; 

lambda = 10;   
epsilon = 0.01;  
tau = 0.001;    

writeout = false;



% experiment_description = 'fixed impulsive noise, rank deficient, uniform, medium conditioned';
experiment_description = 'fixed impulsive noise, rank deficient, uniform, medium conditioned';

txt_output_path = './figures_txt/experiment_II_complex';

data = experiment(m, n, sp, real, lambda, epsilon, tau, maxiter, num_repeats, iter_save, rowsamp, colsamp, ...
           writeout, method_ids, txt_output_path, experiment_description);

















