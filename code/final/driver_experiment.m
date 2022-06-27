% driver for experiment (i) 
% Maximilian Winkler, maximilian.winkler@tu-bs.de
% Lionel Ngoupeyou Tondji, l.ngoupeyou-tondji.tu-bs.de

dir_to_figures = '';
fig_folder_name = '';

method_ids = {'rek', 'srk', 'grek_ad'}; 

m = 5;  % 500 
n = 3;  % 1000
sp = 2;   % 25

real = false;

num_repeats = 2; 
maxiter = 1e4; % iterations

number_data_points = 500;
iter_save = floor(maxiter/number_data_points);  % each such number of iterations, a data point is added in the error plot

rowsamp = 'uniform';
colsamp = 'uniform';

lambda = 5;  
epsilon = 0.01;
tau = 0.001;

writeout = false; 

%experiment_description = 'rank-deficient, medium noise in R(A) complement';

experiment_description = 'rank-deficient, large noise in R(A) complement';


data = experiment(n, m, sp, real, lambda, epsilon, tau, maxiter, num_repeats, iter_save, rowsamp, colsamp,...
                  writeout, method_ids, experiment_description);    
                                         

save('output/data.mat', 'data', '-mat');

                              
                              
                              
                              