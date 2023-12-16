function data = experiment(m, n, sp, real, lambda, epsilon, tau, maxiter, num_repeats, iter_save, ...
                          rowsamp, colsamp, writeout, ...
                          method_ids, txt_output_path, experiment_description)
                          
% Experiment to compare extended Kaczmarz, sparse Kaczmarz and sparse extended Kaczmarz (all randomized)
% We use 
%  g^*(x) = 1/2 ||x||_2^2 or 
%  g^*(x) = r_epsilon(x) + tau/2 ||x||_2^2 
%  with Huber function r_epsilon = ^{epsilon} ||.||_1 (Moreau-envelope of 1-norm) 

% The parameters are:
%  m: number of rows
%  n: number of columns
%  sp: number of nonzeros in exact sparse solution 
%  real: boolean real or complex data A,b, and exact solution x
%  lambda: sparsity parameter in f(x) = 1/2 ||x||_2^2 + lambda ||x||_1
%  epsilon, tau: parameters in g^*(x) = r_epsilon(x) + tau/2 ||x||_2^2, see
%  above for def. of r_epsilon
%  maxiter: maximum number of Kaczmarz-steps 
%  num_repeats: number of repetitions (samples)
%  iter_save: each such number of iterations, a data point is added in the error plot
%  rowsamp (colsamp): Method to sample the rows (columns). Can be
%    'rownorms_squared': probablity proportional to the squares of row-norms
%    ' uniform': Uniform sampling of rows 
%    ' random_probabilities': Sample a vector of probabilities uniformly at
%        random
%  writeout: boolean (if true generate subplot as in paper and save it with
%  matlab2tikz, also save most important data in a data struct,
%  if false show single plots with titles and do not save them)
%  method_ids: cell array with method identifiers. Possible values:
%  'rk': Kaczmarz, 'rek' Extended Kaczmarz, 'srk': Sparse Randomized Kaczmarz, 
%  'gerk_{ad,bd}': (Generalized) Sparse Randomized Extended Kaczmarz,
%  experiment_description: cell array with identifiers for experiments, see
%  'set_up_instance' script

% Some standard parameters for an experiment:
% n = 200;        % Number of columns
% m = 500; % Number of rows
% sp = round(n/20); % Number of nonzeros in solution
% real = true;
% lambda = 5;
% epsilon = 1e-2;
% tau = 1e-3;
% maxiter = 2e5; % Number of iterations
% num_repeats = 50; % Number of repeats over random instances
% iter_save = floor(maxiter/500);
% rowsamp = 'uniform';
% colsamp = 'uniform';
% writeout = false;
% method_ids = {'rek', 'srk', 'gerk_ad'}; 
% experiment_description = 'rank-deficient, large noise in R(A) complement';




addpath('/Users/maximilianwinkler/Documents/Braunschweig/Forschungsthemen/Stochastic_splitting_methods/Kaczmarz method/Sparse Kaczmarz/ExtendedSparseKaczmarz_octave')
addpath('./matlab2tikz')

disp(experiment_description)

fprintf('Experiment A, m=%d, n=%d, sp=%d, maxiter=%d, num_repeats=%d,\n',m,n,sp,maxiter,num_repeats)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper functions for GERK methods
% All used GERK methods are represented in a class and listed below

soft = @(x) max(abs(x)-lambda,0).*sign(x); % Soft thresholding 
grad_Huber = @(x) (1./ max(epsilon,abs(x))+tau).* x; 

gerk_dict = containers.Map;

% GERK-(a,d) method (sparse in x, L2-proj of b onto R(A))
grad_fstar = soft;
L_fstar = 1;
grad_gstar = @(x) x;
L_gstar = 1;
gerk_ad = GERK_method('gerk_ad', grad_fstar, L_fstar, grad_gstar, L_gstar);
gerk_dict('gerk_ad') = gerk_ad;

% GERK-(b,d) method (sparse in x, approximated L1-proj of b onto R(A))
grad_fstar = soft;
L_fstar = 1;
grad_gstar = grad_Huber;
L_gstar = 1/epsilon + tau;
gerk_bd = GERK_method('gerk_bd', grad_fstar, L_fstar, grad_gstar, L_gstar);
gerk_dict('gerk_bd') = gerk_bd;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize storage arrays for results
num_methods = length(method_ids); 
[res, lsres, grad_Huber_functional, err_to_sparse, nonzero_entries]...
                     = deal( zeros(floor(maxiter/iter_save), num_repeats, num_methods) );
x_last_repeat = zeros(n,num_methods);
%[resAbz_list, resATz_list] = deal( zeros(floor(maxiter/iter_save), num_repeats, num_methods) );
%[resAbz_mean_list, resATz_mean_list] = deal( zeros(floor(maxiter/iter_save), num_repeats, num_methods) );
time = zeros(num_repeats, num_methods);
min_rel_cond_A = inf;
max_rel_cond_A = -inf;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up random instances and try methods 

fprintf('experiment ')

for repeats = 1:num_repeats 

    fprintf('# %d, ',repeats)
    
    rand('state',repeats);
    randn('state',repeats);

    problem_data = set_up_instance(m,n,sp,real,experiment_description);
    [A,rel_cond_A,b,bhat,xhat] = deal( problem_data.A, problem_data.rel_cond_A, ...
                                              problem_data.b, problem_data.bhat, problem_data.xhat);   
    min_rel_cond_A = min(rel_cond_A, min_rel_cond_A); 
    max_rel_cond_A = max(rel_cond_A, max_rel_cond_A);
  
    
    norma = sum(abs(A).^2,2);      % squared norms of rows of A
    normc = sum(abs(A).^2,1);  % squared norms of cols of A
    norm_b = norm(b);
    norm_bhat = norm(bhat);
    
    
    % define random functions which do the sampling from the row and column indices 
    
    switch lower(rowsamp)
        case {'rownorms squared'}
            p = norma./sum(norma); 
            P = cumsum(p);
            samp_row = @(k) nnz(rand>P)+1;
        case {'uniform'}
            samp_row = @(k) randi(m,1);           
        case {'random_probabilies'}
            p = rand(m,1); p=p/sum(p); 
            P = cumsum(p);
            samp_row = @(k) nnz(rand>P)+1;
    end     
    
    switch lower(colsamp)
        case {'colnorms squared'}
            pcol = normc./sum(normc);
            Pcol = cumsum(pcol);
            samp_col = @(k) nnz(rand>Pcol)+1;
        case {'uniform'}
            samp_col = @(k) randi(n,1);           
        case {'random_probabilies'}
            pcol = rand(n,1); pcol=pcol/sum(pcol);
            Pcol = cumsum(pcol);
            samp_col = @(k) nnz(rand>Pcol)+1;
    end 
    

    
    
    for method_counter = 1:length(method_ids)
      
        method_id = method_ids{method_counter};
        iter_save_counter = 1;
        tol_nonzero_entries = 1e-5;

        % Initialization
        
        x = zeros(n,1);  
        xstar = zeros(n,1);  
        zstar = b;
        if strcmp(method_id, 'gerk_bd')
           z = gerk_bd.grad_gstar(zstar);
        else
            z = zstar;
        end
        
        
        
        tic;
        
        for iter = 1:maxiter

              %%%%%%%%%%%%
              % 1.: sample indices

              r = samp_row(iter); % tmp - nur zum debuggen
              a = A(r,:);
              
              s = samp_col(iter); % tmp
              c = A(:,s);
              


              %%%%%%%%%%%% 
              % 2.: monitor errors
              
              if mod(iter, iter_save) == 0  
                  
                res(iter_save_counter, repeats, method_counter) = norm(A*x-bhat)/norm_bhat;
                lsres(iter_save_counter, repeats, method_counter) = norm(A'*(A*x-b))/norm_b;
                grad_Huber_functional(iter_save_counter, repeats, method_counter) = norm(A'*gerk_bd.grad_gstar(b-A*x))/norm_b;
                err_to_sparse(iter_save_counter, repeats, method_counter) = norm(x-xhat)/norm(xhat);
                %resAbz_list(iter_save_counter, repeats, method_counter) = norm(A*x-b+zstar)/norm_bhat;
                %resATz_list(iter_save_counter, repeats, method_counter) = norm(A.'*z)/norm_bhat;                
                nonzero_entries(iter_save_counter, repeats, method_counter) = nnz(abs(x) > tol_nonzero_entries);   
                
                iter_save_counter = iter_save_counter + 1;
                
              end
              
              
              %%%%%%%%%%%%
              % 3. perform update 
              
              switch method_id

                % Kaczmarz 
                case 'rk'
                  x = x - (a*x-b(r))/norma(r) *a';

                % Sparse Kaczmarz 
                case 'srk'
                  xstar = xstar -(a*x-b(r))/norma(r)* a';
                  x = soft(xstar);
                  
                % Extended Kaczmarz  
                case 'rek'
                  zstar = zstar - (c'*zstar)/(normc(s)) *c;
                  x = x - (a*x-b(r)+zstar(r))/norma(r) *a';

                % GERK methods
                case {'gerk_ad', 'gerk_bd'}
                  gerk_method = gerk_dict(method_id);
                  zstar = zstar - (c'*z) / (gerk_method.L_gstar * normc(s)) * c;
                  z = gerk_method.grad_gstar(zstar);
                  xstar = xstar -(a*x-b(r)+zstar(r))/(gerk_method.L_fstar * norma(r)) *a';
                  x = gerk_method.grad_fstar(xstar);  

              end 

        end % end for loop over iteration

        % after all iterations
        
        time(repeats, method_counter) = toc;        
        if repeats == num_repeats 
            x_last_repeat(:,method_counter) = x;
        end
        
    end  % end for loop over method

end  % end for loop over repeats








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Settings for plot   

    lightgray =   [0.8 0.8 0.8];
    mediumgray =  [0.6 0.6 0.6];
    lightred =    [1 0.9 0.9];
    mediumred =   [1 0.6 0.6];
    lightgreen =  [0.9 1 0.9];
    mediumgreen = [0.6 1 0.6];
    lightblue =   [0.9 0.9 1];
    mediumblue =  [0.6 0.6 1];
    lightcyan =   [0.9 1 1];
    mediumcyan =  [0.6 1 1];

    linecolor_dict = containers.Map();
    linecolor_dict('rk') = 'c';
    linecolor_dict('srk') = 'b';
    linecolor_dict('rek') = 'k';
    linecolor_dict('gerk_ad') = 'r';
    linecolor_dict('gerk_bd') = 'g';

    minmaxcolor_dict = containers.Map();
    minmaxcolor_dict('rk') = lightcyan;
    minmaxcolor_dict('srk') = lightblue;
    minmaxcolor_dict('rek') = lightgray;
    minmaxcolor_dict('gerk_ad') = lightred;
    minmaxcolor_dict('gerk_bd') = lightgreen;

    quantcolor_dict = containers.Map();
    quantcolor_dict('rk') = mediumcyan;
    quantcolor_dict('srk') = mediumblue;
    quantcolor_dict('rek') = mediumgray;
    quantcolor_dict('gerk_ad') = mediumred; 
    quantcolor_dict('gerk_bd') = mediumgreen; 

    displayname_dict = containers.Map();
    displayname_dict('rk') = 'Kaczmarz';
    displayname_dict('srk') = 'SRK';
    displayname_dict('rek') = 'Extended Kaczmarz';
    displayname_dict('gerk_ad') = 'GERK-(a,d)'; 
    displayname_dict('gerk_bd') = 'GERK-(b,d)'; 

    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up error arrays for plot

    [min_res, max_res, median_res, quant25_res, quant75_res] = compute_minmax_median_quantiles(res);
    [min_lsres, max_lsres, median_lsres, quant25_lsres, quant75_lsres] = compute_minmax_median_quantiles(lsres);
    [min_err_to_sparse, max_err_to_sparse, median_err_to_sparse, quant25_err_to_sparse, quant75_err_to_sparse]...
         = compute_minmax_median_quantiles(err_to_sparse);
    [min_nonzero_entries, max_nonzero_entries, median_nonzero_entries, ~, ~]...
         = compute_minmax_median_quantiles(nonzero_entries);
    [min_grad_Huber_functional, max_grad_Huber_functional, median_grad_Huber_functional, quant25_grad_Huber_functional,...
     quant75_grad_Huber_functional]...
         = compute_minmax_median_quantiles(grad_Huber_functional);  
    
    
    % set zero entries to eps for not getting -inf in log plot
    min_res = set_zero_entries_to_eps(min_res); 
    min_lsres = set_zero_entries_to_eps(min_lsres); 
    min_grad_Huber_functional = set_zero_entries_to_eps(min_grad_Huber_functional);
    min_err_to_sparse = set_zero_entries_to_eps(min_err_to_sparse);

   
    
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % make plots
    
    % plot no. 1: errors
    
    
    if any(strcmp(method_ids, 'gerk_bd'))
        num_subplots = 4;
    else
        num_subplots = 3;
    end


    figure

    % ||Ax- b_{exact}|| / ||b_{exact}|| 
    sgtitle(sprintf('Errors, m = %d, n = %d, s = %d, repeats = %d',m,n,sp,num_repeats));
    h1 = subplot(1,num_subplots,1);
    choose_logy = true;
    [plot_array, ~, ~] = plot_minmax_median_quantiles('-',min_res,max_res,median_res,quant25_res,quant75_res,choose_logy,method_ids,iter_save,maxiter,minmaxcolor_dict,quantcolor_dict,linecolor_dict,displayname_dict);
    if ~ writeout
        legend(plot_array, 'location', 'northwest', 'FontSize', 5);
    end

    xlabel('$k$','Interpreter','latex')
    ylabel('$$\|Ax-\hat b\|/\|\hat b\|$$','Interpreter','latex')





    % ||A^T(b-Ax)||
    h2 = subplot(1,num_subplots,2);
    choose_logy = true;
    plot_minmax_median_quantiles('-',min_lsres,max_lsres,median_lsres,...
                                                quant25_lsres,quant75_lsres,choose_logy, ...
                                                method_ids,iter_save,maxiter,minmaxcolor_dict,quantcolor_dict,...
                                                linecolor_dict,displayname_dict);
    xlabel('$k$','Interpreter','latex')
    ylabel('$$\|A^T(b-Ax)\|/\|b\|$$','Interpreter','latex')     



    % ||A^T\nabla g^*(b-Ax)||
    if any(strcmp(method_ids, 'gerk_bd'))            
        h3 = subplot(1,num_subplots,3);
        choose_logy = true;
        plot_minmax_median_quantiles('-',min_grad_Huber_functional,max_grad_Huber_functional,median_grad_Huber_functional,...
                                                    quant25_grad_Huber_functional,quant75_grad_Huber_functional,choose_logy, ...
                                                    method_ids,iter_save,maxiter,minmaxcolor_dict,quantcolor_dict,...
                                                    linecolor_dict,displayname_dict);
        xlabel('$k$','Interpreter','latex')
        ylabel('$$\|A^T \nabla g^*(b-Ax)\|/\|b\|$$','Interpreter','latex')        
    end




    % || x-\hat x||   
    h4 = subplot(1,num_subplots,num_subplots);
    choose_logy = true;
    plot_minmax_median_quantiles('-',min_err_to_sparse,max_err_to_sparse,median_err_to_sparse,...
                                                    quant25_err_to_sparse,quant75_err_to_sparse, choose_logy,...
                                                    method_ids,iter_save,maxiter,minmaxcolor_dict,...
                                                    quantcolor_dict,linecolor_dict,displayname_dict);

    xlabel('$k$','Interpreter','latex')
    ylabel('$$\|x-\hat x\|/\|\hat x\|$$','Interpreter','latex') 

    



    % add horizontal space between subplots
    space = 0.05;
    pos_h1 = get(h1,'Position');  % inner box of first subplot (left, bottom, width, height)
    pos_h1(1) = pos_h1(1) - space;

    set(h1, 'Position', pos_h1);
    drawnow;


    if any(strcmp(method_ids, 'gerk_bd'))

        pos = get(h2, 'Position');        
        set(h2, 'Position', [pos(1)-0.3*space, pos_h1(2:end)]);
        drawnow;

        pos = get(h3, 'Position');        
        set(h3, 'Position', [pos(1)+0.4*space, pos_h1(2:end)]);
        drawnow;

        pos = get(h4, 'Position');
        set(h4, 'Position', [pos(1)+1.2*space, pos_h1(2:end)]);    
        drawnow;

    else
        
        pos = get(h2, 'Position');        
        set(h2, 'Position', [pos(1)-0.3*space, pos_h1(2:end)]);
        drawnow;

        pos = get(h4, 'Position');        
        set(h4, 'Position', [pos(1)+0.4*space, pos_h1(2:end)]);
        drawnow;
     
    end



    % save first row of figures with matlab2tikz
    matlab2tikz('width','\figurewidth',...
    'output/err_over_iter.tex');        


    % save first row of figures in txt file (errors over iterations)
%{
    num_iter_array = (1:iter_save:maxiter)';

    % (1,1)
    writeout_data_over_array_on_xaxis(...
    [txt_output_path '/res_over_iter.txt'],...
    {'k', 'min_rek', 'min_srk', 'min_gerk_ad', 'min_gerk_bd', ...
          'max_rek', 'max_srk', 'max_gerk_ad', 'max_gerk_bd', ...
          'median_rek', 'median_srk', 'median_gerk_ad', 'median_gerk_bd',...
          'quant25_rek', 'quant25_srk', 'quant25_gerk_ad', 'quant25_gerk_bd',...
          'quant75_rek', 'quant75_srk', 'quant75_gerk_ad', 'quant75_gerk_bd'}, ...
    num_iter_array, ...
    [min_res(:,1), min_res(:,2), min_res(:,3), ...
     max_res(:,1), max_res(:,2), max_res(:,3), ...
     median_res(:,1), median_res(:,2), median_res(:,3), ...
     quant25_res(:,1), quant25_res(:,2), quant25_res(:,3), ...
     quant75_res(:,1), quant75_res(:,2), quant75_res(:,3)]); 

    % (1,2)
    writeout_data_over_array_on_xaxis(...
    [txt_output_path '/lsres_over_iter.txt'],...
    {'k', 'min_rek', 'min_srk', 'min_gerk_ad', 'min_gerk_bd', ...
          'max_rek', 'max_srk', 'max_gerk_ad', 'max_gerk_bd', ...
          'median_rek', 'median_srk', 'median_gerk_ad', 'median_gerk_bd', ...
          'quant25_rek', 'quant25_srk', 'quant25_gerk_ad', 'quant25_gerk_bd', ...
          'quant75_rek', 'quant75_srk', 'quant75_gerk_ad', 'quant75_gerk_bd'}, ...
    num_iter_array, ...
    [min_lsres(:,1), min_lsres(:,2), min_lsres(:,3), ...
     max_lsres(:,1), max_lsres(:,2), max_lsres(:,3), ...
     median_lsres(:,1), median_lsres(:,2), median_lsres(:,3), ...
     quant25_lsres(:,1), quant25_lsres(:,2), quant25_lsres(:,3), ...
     quant75_lsres(:,1), quant75_lsres(:,2), quant75_lsres(:,3)]); 

    % (1,3)
    writeout_data_over_array_on_xaxis(...
    [txt_output_path '/grad_Huber_functional.txt'],...
    {'k', 'min_rek', 'min_srk', 'min_gerk_ad', 'min_gerk_bd', ...
          'max_rek', 'max_srk', 'max_gerk_ad', 'max_gerk_bd', ...
          'median_rek', 'median_srk', 'median_gerk_ad', 'median_gerk_bd',...
          'quant25_rek', 'quant25_srk', 'quant25_gerk_ad', 'quant25_gerk_bd',...
          'quant75_rek', 'quant75_srk', 'quant75_gerk_ad', 'quant75_gerk_bd'}, ...
    num_iter_array, ...
    [min_grad_Huber_functional(:,1), min_grad_Huber_functional(:,2), min_grad_Huber_functional(:,3), ...
    max_grad_Huber_functional(:,1), max_grad_Huber_functional(:,2), max_grad_Huber_functional(:,3), ...
    median_grad_Huber_functional(:,1), median_grad_Huber_functional(:,2), median_grad_Huber_functional(:,3), ...
    quant25_grad_Huber_functional(:,1), quant25_grad_Huber_functional(:,2), quant25_grad_Huber_functional(:,3), ...
    quant75_grad_Huber_functional(:,1), quant75_grad_Huber_functional(:,2), quant75_grad_Huber_functional(:,3)]); 

    % (1,4)
    writeout_data_over_array_on_xaxis(...
    [txt_output_path '/err_to_sparse_over_iter.txt'],...
    {'k', 'min_rek', 'min_srk', 'min_gerk_ad', 'min_gerk_bd', ...
          'max_rek', 'max_srk', 'max_gerk_ad', 'max_gerk_bd', ...
          'median_rek', 'median_srk', 'median_gerk_ad', 'median_gerk_bd',...
          'quant25_rek', 'quant25_srk', 'quant25_gerk_ad', 'quant25_gerk_bd',...
          'quant75_rek', 'quant75_srk', 'quant75_gerk_ad', 'quant75_gerk_bd'}, ...
    num_iter_array, ...
    [min_err_to_sparse(:,1), min_err_to_sparse(:,2), min_err_to_sparse(:,3), ...
     max_err_to_sparse(:,1), max_err_to_sparse(:,2), max_err_to_sparse(:,3), ...
     median_err_to_sparse(:,1), median_err_to_sparse(:,2), median_err_to_sparse(:,3), ...
     quant25_err_to_sparse(:,1), quant25_err_to_sparse(:,2), quant25_err_to_sparse(:,3), ...
     quant75_err_to_sparse(:,1), quant75_err_to_sparse(:,2), quant75_err_to_sparse(:,3)]); 
%}

    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot no. 2: stem plots of corrupted right-hand side 
    % and last iterates x together with exact sparse solution xhat
 
    figure

    num_subplots = 1 + length(method_ids); 

    % b_{exact} vs. b       

    subplot(1,num_subplots,1);

    if real 
        plot_array = plot(1:m,bhat,1:m,b,':'); 
    else
        plot_array = plot(1:m,abs(bhat),1:m,abs(b),':'); 
    end

    axis square

    % remove legend and xlabels
    leg = legend(plot_array, {'b','noisy b'}, 'location', 'northwest', 'FontSize', 7);
    if writeout
        set(leg, 'visible', 'off');
        set(gca,'XTick',[])
    end


    plot_nr = 2;

    for method_counter = 1:length(method_ids) 

        subplot(1,num_subplots,plot_nr);

        % only plot last iterate
        x = x_last_repeat(:,method_counter);
        
        if ~real  % plot absolute values
            xhat = abs(xhat);
            x = abs(x);
        end

        stem(xhat,'color','blue')
        hold on
        stem(x,'color','red')
        axis square
        
       leg = legend({'x_{hat}','x'}, 'location', 'northwest', 'FontSize', 7);
       if writeout
           set(leg, 'visible', 'off');
           set(gca,'XTick',[])
       end

        plot_nr = plot_nr + 1;

    end

    % save figure with matlab2tikz
    matlab2tikz('width','\figurewidth',...
    'output/x_components.tex');




    % save second row of figures in txt file (stem plots of solutions)
%{
    xhat_plot = xhat;
    x_plot = x_last_repeat;
    if ~real
        xhat_plot = abs(xhat_plot);
        x_plot = abs(x_plot);
        bhat = abs(bhat);
        b = abs(b);
    end

    writeout_data_over_array_on_xaxis(...
    [txt_output_path '/b_and_noisy_b.txt'],...
    {'i', 'bhat', 'b'}, ...
    1:length(b), ...
    [bhat, b]); 

    writeout_data_over_array_on_xaxis(...
    [txt_output_path '/x_true_and_approx.txt'],...
    {'i', 'xhat', 'x_rek', 'x_srk', 'x_gerk_ad', 'x_gerk_bd'}, ...
    1:length(xhat_plot), ...
    [xhat_plot, x_plot(:,1), x_plot(:,2), x_plot(:,3)]); 
%}




    % Display computation time and number of nonzero entries
    
    fprintf('\n')
    for i = 1:length(method_ids)
        disp(['Entries with absolute value > ' num2str(tol_nonzero_entries) ', ' method_ids{i} ' method (min, median,max): ' num2str(min_nonzero_entries(end,i)) ', ' num2str(median_nonzero_entries(end,i)) ', ' num2str(max_nonzero_entries(end,i)) ])
    end

    min_time = min(time,[],1);
    max_time = max(time,[],1);
    median_time = median(time,1);
    for i = 1:length(method_ids)
      disp(['Computation time for ' num2str(maxiter) ' iterations, ' method_ids{i} ' method(min,median,max): ' num2str(min_time(i)) ', ' num2str(median_time(i)) ', ' num2str(max_time(i))])
    end





    % save data

    data = struct('methods', method_ids, 'min_rel_cond_A', min_rel_cond_A, 'max_rel_cond_A', max_rel_cond_A, ...
                  'x_last_repeat', x_last_repeat, ...
                  'max_res', max_res, 'min_res', min_res, 'median_res', median_res, ...
                  'quant25_res', quant25_res, 'quant75_res', quant75_res,...
                  'max_lsres', max_lsres, 'min_lsres', min_lsres, 'median_lsres', median_lsres, ...
                  'quant25_lsres', quant25_lsres, 'quant75_lsres', quant75_lsres,...
                  'max_err_to_sparse', max_err_to_sparse, 'min_err_to_sparse', min_err_to_sparse, 'median_err_to_sparse', median_err_to_sparse, ...
                  'quant25_err_to_sparse', quant25_err_to_sparse, 'quant75_err_to_sparse', quant75_err_to_sparse,...
                  'max_nonzero_entries', max_nonzero_entries, 'min_nonzero_entries', min_nonzero_entries, ...
                  'median_nonzero_entries', median_nonzero_entries);

    if writeout
        save('output/data.mat', 'data', '-mat');
    end
              
              
         
              
              
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
    % helper functions for visualization
    

    % for logarithmic plot
    function arr = set_zero_entries_to_eps(arr)

        arr(arr<eps) = eps;

    end




    function [mins, maxs, medians, quant25s, quant75s] = compute_minmax_median_quantiles(arr)
        for k = 1:size(arr,3)  % loop over all methods

          if size(arr(:,:,k),2) == 1   % if only one repeat, nothing to compute 
            [mins,maxs,medians,quant25s,quant75s] = deal(arr(:,1,:));
          else
            mins(:,k) = min(arr(:,:,k), [], 2);
            maxs(:,k) = max(arr(:,:,k), [], 2);
            medians(:,k) = median(arr(:,:,k),2);
            quant25s(:,k) = quantile(arr(:,:,k),0.25,2); 
            quant75s(:,k) = quantile(arr(:,:,k),0.75,2);
          end
        end
    end







end




