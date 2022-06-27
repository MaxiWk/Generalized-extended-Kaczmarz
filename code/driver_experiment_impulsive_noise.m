
method_ids = {'rek','srk','gerk_ad','gerk_bd'}; 

m = 3;   % 1000
n = 5;  % 500
sp = 2;  % min(m,n)/20; 

real = false; 

num_repeats = 2; % 50
maxiter = 2e5; %2e5 

number_data_points = 500;
iter_save = floor(maxiter/number_data_points);  % each such number of iterations, a data point is added in the error plot

rowsamp = 'uniform';
colsamp = 'uniform'; 

lambda = 5;   
epsilon = 0.01;  
tau = 0.001;    

writeout = true;



%experiment_description = 'impulsive noise, rank deficient, uniform, well conditioned';
%experiment_description = 'impulsive noise, rank deficient, uniform, bad conditioned';
%experiment_description = 'impulsive noise, rank deficient, uniform, medium conditioned';
experiment_description = 'fixed impulsive noise, rank deficient, uniform, medium conditioned';

data = experiment(m, n, sp, real, lambda, epsilon, tau, maxiter, num_repeats, iter_save, rowsamp, colsamp, ...
           writeout, method_ids, experiment_description);

















