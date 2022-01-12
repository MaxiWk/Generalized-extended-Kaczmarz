function data = experiment(n,m,sp,real_setting,lambda,T,L_gstar,maxiter,num_repeats,iter_save,rowsamp,colsamp,writeout,dir_to_figures,fig_folder_name,disp_instance,savestep,stopcrit_sample_pars,method_array,experiment_description)
% Experiment to ompare usual extended Kaczmarz, sparse Kaczmarz and sparse extended Kaczmarz (all randomized)
% The parameters are:
%  n: number of columns
%  m: number of rows
%  sp: number of nozeros in solution 
%  real_setting: boolean, if true: solve Ax=b with real data A,b, if false:
%  complex data A,b
%  T: gradient g*, can also be a cell array of different gradients g* -> then names 'grek_1', 'grek_2',.. 
%  Needs to match also to L_gstar_1, L_gstar_2, ...
%  We require length(T) == length(L_gstar) < 3 if cell array
%  We use 
%  g^*(x) = 1/2 ||x||_2^2 or 
%  g^*(x) = r_epsilon(x) + tau/2 ||x||_2^2 with Huber function r_epsilon =
%  ^{epsilon} ||.||_1 (Moreau-envelope of 1-norm)
%  Here always gradient f^* = soft shrinkage with parameter lambda
%  L_gstar: Lipschitz constant of T (i.e. of gradient g^*)
%  maxiter: maximum number of Kaczmarz-steps (projections, not sweeps)
%  num_repeats: number of repetitions (samples)
%  iter_save: each such number of iterations, a data point is added in the error plot
%  rowsamp: Method to sample the rows. Can be
%    'rownorms_squared': probablity proportional to the squares of row-norms
%    ' uniform': Uniform sampling of rows 
%    ' random_probabilities': Sample a vector of probabilities uniformly at
%        random
%  colsamp: Analogously for column sampling in case of extended Kaczmarz
%  writeout: boolean (if true generate subplot as in paper and save it with matlab2tikz, if false show single plots with titles and do not save them)
%  dir_to_figures: Directory to folder where to save figures with matlab2tikz
%  in case of writeout==true
%  fig_folder_name: in case of writeout==true: is added on dir_to_figures and in the name of the error
%  disp_instance: boolean (if true disp b_exact, the right-hand side without noise, and xhat, the exact sparse solution)
%  savestep: Only each "savesteps"th data point will be stored
%  stopcrit_sample_pars: struct with tolerances for stopping criterion, see set_up_instance
%  script for examples
%  method_array: cell array with method identifiers. Possible values:
%  'rk': Kaczmarz, 'srk': Sparse Randomized Kaczmarz, 'esrk': Exact step Sparse Kaczmarz, 
%  'grek_{1,2,3}': (Generalized) Sparse Randomized Extended Kaczmarz, 'egrek': Exact step Sparse
%  Randomized Extended Kaczmarz
%  experiment_description: cell array with identifiers for experiments, see
%  set_up_instance script

% Some standard parameters for an experiment:
% n = 200;        % Number of columns
% m = round(5*n); % Number of rows
% sp = round(n/8); % Number of nonzeros in solution
% maxiter = 40*m; % Number of iterations
% rowsamp = 'rownorms squared';
% num_repeats = 60;% Number of repeats over random instances
% writeout = false;



addpath('/Users/maximilianwinkler/Documents/Braunschweig/Forschungsthemen/Stochastic_splitting_methods/Kaczmarz method/Sparse Kaczmarz/ExtendedSparseKaczmarz_octave')
addpath('./matlab2tikz')

disp(experiment_description)

fprintf('Experiment A, m=%d, n=%d, sp=%d, maxiter=%d, num_repeats=%d,\n',m,n,sp,maxiter,num_repeats)



% 'unpack' T and L_gstar list 
for i = 1:length(T)
    if i>0
        T_1 = T{1};
        L_gstar_1 = L_gstar(1);
    end
    if i>1
        T_2 = T{2};
        L_gstar_2 = L_gstar(2);
    end
    if i>2
        T_3 = T{3};
        L_gstar_3 = L_gstar(3);
    end
end



% Initialize storage arrays for results
num_methods = length(method_array); 
[res, lsres, grad_zfunctional, err_to_sparse, err_to_moorepi, nonzero_entries]...
                     = deal( zeros(floor(maxiter/iter_save), num_repeats, num_methods) );
err_to_sparse_stopped = zeros(num_repeats, num_methods);   
err_to_sparse_end_div_stopped = zeros(num_repeats, num_methods); 
x_last_test_instance = zeros(n,num_methods);
[length_resAbz_sampled, length_resATz_sampled] = deal(stopcrit_sample_pars.length_resAbz_sampled,...
                                                             stopcrit_sample_pars.length_resATz_sampled);
resAbz_sampled = zeros(1, stopcrit_sample_pars.length_resAbz_sampled);
resATz_sampled = zeros(1, stopcrit_sample_pars.length_resATz_sampled);
[resAbz_mean, resATz_mean] = deal( zeros(1,maxiter) );
iterstop_list = zeros(num_repeats,num_methods);
xstop_list = zeros(n, num_repeats, num_methods);
[resAbz_list, resATz_list] = deal( zeros(floor(maxiter/iter_save), num_repeats, num_methods) );
[resAbz_mean_list, resATz_mean_list] = deal( zeros(floor(maxiter/iter_save), num_repeats, num_methods) );
time = zeros(num_repeats, num_methods);

min_rel_cond_A = inf;
max_rel_cond_A = -inf;





% Loop over number of repeats 
fprintf('experiment ')
for repeats = 1:num_repeats 

    fprintf('# %d, ',repeats)
    
    rand('state',repeats);
    randn('state',repeats);

    problem_data = set_up_instance(m,n,sp,real_setting,experiment_description);
    [A,rel_cond_A,b,b_exact,xhat] = deal( problem_data.A, problem_data.rel_cond_A, ...
                                              problem_data.b, problem_data.b_exact, problem_data.xhat);
    [tol_resAbz, tol_resATz] = deal( problem_data.tol_resAbz, problem_data.tol_resATz );
    %rankA = rank(A);
    
    min_rel_cond_A = min(rel_cond_A, min_rel_cond_A); 
    max_rel_cond_A = max(rel_cond_A, max_rel_cond_A);
    
    if disp_instance 
      disp('')
      disp('xhat (ground truth): ')
      disp(xhat)
      disp('b (ground truth)')
      disp(b_exact)
    end
  
    xmoorepi = pinv(A)*b;
    %proj_b = A*xmoorepi;
    
    
    
    S = @(z) max(abs(z)-lambda,0).*sign(z); % Soft thresholding    

    
    norma = sum(abs(A).^2,2);      % squared norms of rows of A
    norma_col = sum(abs(A).^2,1);  % squared norms of cols of A
    norm_b_exact = norm(b_exact);
    
    
    % define functions which sample from the row and column indices 
    
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
    end % end row sampling switch
    
    switch lower(colsamp)
        case {'colnorms squared'}
            pcol = norma_col./sum(norma_col);
            Pcol = cumsum(pcol);
            samp_col = @(k) nnz(rand>Pcol)+1;
        case {'uniform'}
            samp_col = @(k) randi(n,1);           
        case {'random_probabilies'}
            pcol = rand(n,1); pcol=pcol/sum(pcol);
            Pcol = cumsum(pcol);
            samp_col = @(k) nnz(rand>Pcol)+1;
    end % end row sampling switch
    
    last_repeat = (repeats == num_repeats);
   
    
    % tmp (nur zum debuggen)
    row_ind=randi([1,m],1,maxiter);
    col_ind=randi([1,n],1,maxiter);
    
    
    
    for method_counter = 1:length(method_array)
      
        method = method_array{method_counter};
        iter_save_counter = 1;

        % Initialize method
        
        x = zeros(n,1);  
        % tmp start
        %x0 = randn(size(A,2),1);
        %x = x0;
        %P_Axb_x0 = x - pinv(A)*(A*x-b);
        % tmp end
        xdual = zeros(n,1);  
        zdual = b;
        z = T_1(zdual);
        
        
        % precompute expressions for deterministic algorithms (only for small dimensions)
        
        if any(strcmp(method_array,'det_sek'))
          AAT = A*A.'; 
        end
        
        if any(strcmp(method_array,'lin_breg'))
          ATA = A.'*A;
          ATb = A.'*b;
        end
        
        lin_breg_stepsize = 1/norm(A)^2;
        
        tic;
        

        
        stopped = false;
        
        for iter = 1:maxiter

              %%%%%%%%%%%%
              % 1.: choose index samples

              r = samp_row(iter);
              r = row_ind(iter); % tmp - nur zum debuggen
              a = A(r,:);
              
              s = samp_col(iter); % execute this for all algorithms in order to get the same r indices
              s = col_ind(iter); % tmp
              c = A(:,s);
              
              
              %%%%%%%%%%%% 
              % 2.: check stopping criterion

              resAbz_sampled(iter) = (a*x-b(r)+zdual(r));
              first_iter_with_row_samples = max(1, iter-length_resAbz_sampled);
              actual_length_resAbz_sampled = iter - first_iter_with_row_samples +1;
              resAbz_mean(iter) = norm(resAbz_sampled(first_iter_with_row_samples:iter)) *sqrt(m/max(1,actual_length_resAbz_sampled));
              
              resATz_sampled(iter) = c'*z;
              first_iter_with_col_samples = max(1,iter-length_resATz_sampled);
              actual_length_resATz_sampled = iter - first_iter_with_col_samples +1;
              resATz_mean(iter) = norm(resATz_sampled(first_iter_with_col_samples:iter)) * sqrt(n/max(1,actual_length_resATz_sampled));             

              if( ~stopped ...
                  && iter > stopcrit_sample_pars.min_possible_iter_for_stopping ...
                  && resAbz_mean(iter) < tol_resAbz ...
                  && resATz_mean(iter) < norm(b-zdual)*tol_resATz )

                  iterstop_list(repeats,method_counter) = iter;
                  xstop_list(:,repeats,method_counter) = x;
                  err_to_sparse_stopped(repeats,method_counter) = norm(x-xhat)/norm(xhat);                 
                  stopped = true;                  
                  
                  %break; 
                  %err_to_sparse_end_div_stopped(repeats,method_counter) = 1;

              end


              %%%%%%%%%%%% 
              % 3.: save error quntities 
              
              if mod(iter, iter_save) == 0  
                %disp(num2str(norm(x-xmoorepi)));  %%%%%%%
                res(iter_save_counter, repeats, method_counter) = norm(A*x-b_exact)/norm_b_exact;
                %lsres(iter_save_counter, repeats, method_counter) = norm(A'*(A*x-A*xhat))/norm_b_exact;
                %res_proj(iter_save_counter, repeats, method_counter) = norm(A*x - proj_b)/norm_b_exact;
                lsres(iter_save_counter, repeats, method_counter) = norm(A.'*(A*x-b_exact))/norm_b_exact;
                %lsres_proj(iter_save_counter, repeats, method_counter) = norm(A'*(A*x - proj_b))/norm_ATb;
                grad_zfunctional(iter_save_counter, repeats, method_counter) = norm(A.'*T_1(b_exact-A*x))/norm_b_exact;
                %err_to_moorepi(iter_save_counter, repeats, method_counter) = norm(x-xmoorepi)/norm(xmoorepi);
                err_to_sparse(iter_save_counter, repeats, method_counter) = norm(x-xhat)/norm(xhat);
                resAbz_list(iter_save_counter, repeats, method_counter) = norm(A*x-b+zdual)/norm_b_exact;
                resATz_list(iter_save_counter, repeats, method_counter) = norm(A.'*z)/norm_b_exact;
                resAbz_mean_list(iter_save_counter, repeats, method_counter) = resAbz_mean(iter);
                resATz_mean_list(iter_save_counter, repeats, method_counter) = resATz_mean(iter);
                tol_zero = 1e-5;   % we count entries with abs value > tol_zero
                nonzero_entries(iter_save_counter, repeats, method_counter) = nnz(abs(x) > tol_zero);        
                iter_save_counter = iter_save_counter + 1;
                
              end
              
              
              %%%%%%%%%%%%
              % 4. perform update 
              
              switch method

                % Kaczmarz 
                case 'rk'
                  x = x - (a*x-b(r))/norma(r) *a';

                % Sparse Kaczmarz 
              case 'srk'
                  xdual = xdual -(a*x-b(r))/norma(r)* a';
                  x = S(xdual);

                % Sparse Kaczmarz with exact stept
                case 'esrk'
                  [x,xdual] = linesearch_shrinkage(x,xdual,a.',b(r),lambda);

                % Extended Kaczmarz  
                case 'rek'
                  zdual = zdual - (c'*zdual)/(norma_col(s)) *c;
                  x = x - (a*x-b(r)+zdual(r))/norma(r) *a';

                % Sparse Extended Kaczmarz 
                case 'grek_1'
                  zdual = zdual - (c'*z)/(L_gstar_1*norma_col(s)) *c;
                  z = T_1(zdual);
                  xdual = xdual -(a*x-b(r)+zdual(r))/norma(r) *a';
                  x = S(xdual);  

                % Sparse Extended Kaczmarz 
                case 'grek_2'
                  zdual = zdual - (c'*z)/(L_gstar_2*norma_col(s)) *c;   
                  z = T_2(zdual);
                  xdual = xdual -(a*x-b(r)+zdual(r))/norma(r) *a';
                  x = S(xdual);  
                  
                % Sparse Extended Kaczmarz 
                case 'grek_3'
                  zdual = zdual - (c'*z)/(L_gstar_3*norma_col(s)) *c;
                  z = T_3(zdual);
                  xdual = xdual -(a*x-b(r)+zdual(r))/norma(r) *a';
                  x = S(xdual);                    
                  
                % Sparse Extended Kaczmarz with exact step
                case 'egrek'
                  zdual = zdual - (c'*z)/(L_gstar*norma_col(s)) *c;
                  z = T_1(zdual);
                  [x,xdual] = linesearch_shrinkage(x,xdual,a.',b(r)-zdual(r),lambda);

                case 'det_sek'
                  zdual = zdual - lin_breg_stepsize*AAT*zdual;
                  z = T(zdual);
                  xdual = xdual - lin_breg_stepsize*A.'*(A*x-b+zdual);
                  x = S(xdual);

              case 'lin_breg'
                  xdual = xdual - lin_breg_stepsize*(ATA*x-ATb);
                  x = S(xdual);

              end 

        end % end for loop over iteration

         

        if ~stopped
            iterstop_list(repeats,method_counter) = maxiter;
            xstop_list(:,repeats,method_counter) = x;
            err_to_sparse_stopped(repeats,method_counter) = err_to_sparse(end,repeats,method_counter);
            err_to_sparse_end_div_stopped(repeats,method_counter) = 1;
        end
        
        % after all iterations
        time(repeats, method_counter) = toc;
        
        if last_repeat 
            x_last_test_instance(:,method_counter) = x;
        end
        
        err_to_sparse_end_div_stopped(repeats,method_counter)...
            = err_to_sparse(end,repeats,method_counter) / err_to_sparse_stopped(repeats,method_counter);
    
    end  % end for loop over method

end  % end for loop over repeats








    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Plot results and export data for residuals

    %close all

    lightgray =   [0.8 0.8 0.8];
    mediumgray =  [0.6 0.6 0.6];
    lightred =    [1 0.9 0.9];
    mediumred =   [1 0.6 0.6];
    lightgreen =  [0.9 1 0.9];
    mediumgreen = [0.6 1 0.6];
    lightblue =   [0.9 0.9 1];
    mediumblue =  [0.6 0.6 1];
    lightyellow = [1 1 0.9];
    mediumyellow =[1 1 0.6];
    lightcyan =   [0.9 1 1];
    mediumcyan =  [0.6 1 1];

    linecolor_dict = containers.Map();
    linecolor_dict('rk') = 'c';
    linecolor_dict('srk') = 'b';
    linecolor_dict('esrk') = 'y';
    linecolor_dict('rek') = 'k';
    linecolor_dict('grek_1') = 'r';
    linecolor_dict('grek_2') = 'g';
    linecolor_dict('grek_3') = 'y';
    linecolor_dict('egrek') = 'g';
    %linecolor_dict('det_sek') = 'c';
    linecolor_dict('lin_breg') = 'k';

    minmaxcolor_dict = containers.Map();
    minmaxcolor_dict('rk') = lightcyan;
    minmaxcolor_dict('srk') = lightblue;
    minmaxcolor_dict('esrk') = lightyellow;
    minmaxcolor_dict('rek') = lightgray;
    minmaxcolor_dict('grek_1') = lightred;
    minmaxcolor_dict('grek_2') = lightgreen;
    minmaxcolor_dict('grek_3') = lightgray;
    minmaxcolor_dict('egrek') = lightgreen;
    %minmaxcolor_dict('det_sek') = lightcyan;
    minmaxcolor_dict('lin_breg') = lightgray; 

    quantcolor_dict = containers.Map();
    quantcolor_dict('rk') = mediumcyan;
    quantcolor_dict('srk') = mediumblue;
    quantcolor_dict('esrk') = mediumyellow;
    quantcolor_dict('rek') = mediumgray;
    quantcolor_dict('grek_1') = mediumred; 
    quantcolor_dict('grek_2') = mediumgreen; 
    quantcolor_dict('grek_3') = mediumgray; 
    quantcolor_dict('egrek') = mediumgreen;
    %quantcolor_dict('det_sek') = mediumcyan; 
    quantcolor_dict('lin_breg') = mediumgray;

    displayname_dict = containers.Map();
    displayname_dict('rk') = 'Kaczmarz';
    displayname_dict('srk') = 'SRK';
    displayname_dict('esrk') = 'ESRK';
    displayname_dict('rek') = 'Extended Kaczmarz';
    displayname_dict('grek_1') = 'GREK 1'; 
    displayname_dict('grek_2') = 'GREK 2'; 
    displayname_dict('grek_3') = 'GREK 3';
    displayname_dict('egrek') = 'GREK with ESRK stepsize';
    displayname_dict('lin_breg') = 'Linearized Bregman method';


    %% Set up plotable arrays (min, max, quantiles over all repeats)


    % if only one repeat comment out from here ------

    [min_res, max_res, median_res, quant25_res, quant75_res] = compute_minmax_median_quantiles(res);
    [min_lsres, max_lsres, median_lsres, quant25_lsres, quant75_lsres] = compute_minmax_median_quantiles(lsres);
    %[min_res_proj, max_res_proj, median_res_proj, quant25_res_proj, quant75_res_proj] = compute_minmax_median_quantiles(res_proj);
    %[min_lsres_proj, max_lsres_proj, median_lsres_proj, quant25_lsres_proj, quant75_lsres_proj] = compute_minmax_median_quantiles(lsres_proj);
    [min_err_to_moorepi, max_err_to_moorepi, median_err_to_moorepi, quant25_err_to_moorepi, quant75_err_to_moorepi]...
         = compute_minmax_median_quantiles(err_to_moorepi);
    [min_err_to_sparse, max_err_to_sparse, median_err_to_sparse, quant25_err_to_sparse, quant75_err_to_sparse]...
         = compute_minmax_median_quantiles(err_to_sparse);
    [min_nonzero_entries, max_nonzero_entries, median_nonzero_entries, quant25_nonzero_entries, quant75_nonzero_entries]...
         = compute_minmax_median_quantiles(nonzero_entries);
    [min_grad_zfunctional, max_grad_zfunctional, median_grad_zfunctional, quant25_grad_zfunctional,...
     quant75_grad_zfunctional]...
         = compute_minmax_median_quantiles(grad_zfunctional);  
    [min_resAbz_list, max_resAbz_list, median_resAbz_list, quant25_resAbz_list, quant75_resAbz_list] = compute_minmax_median_quantiles(resAbz_list);
    [min_resATz_list, max_resATz_list, median_resATz_list, quant25_resATz_list, quant75_resATz_list] = compute_minmax_median_quantiles(resATz_list);
    [min_resAbz_mean_list, max_resAbz_mean_list, median_resAbz_mean_list, quant25_resAbz_mean_list, quant75_resAbz_mean_list] = compute_minmax_median_quantiles(resAbz_mean_list);
    [min_resATz_mean_list, max_resATz_mean_list, median_resATz_mean_list, quant25_resATz_mean_list, quant75_resATz_mean_list] = compute_minmax_median_quantiles(resATz_mean_list);

    
    
    % set zero entries to eps for not getting -inf in log plot
    
    min_res = set_zero_entries_to_eps(min_res); 
    %min_res_proj = set_zero_entries_to_eps(min_res_proj); 
    min_lsres = set_zero_entries_to_eps(min_lsres); 
    min_grad_zfunctional = set_zero_entries_to_eps(min_grad_zfunctional);
    %min_lsres_proj = set_zero_entries_to_eps(min_lsres_proj);    
    min_err_to_moorepi = set_zero_entries_to_eps(min_err_to_moorepi);
    min_err_to_sparse = set_zero_entries_to_eps(min_err_to_sparse);




    
    if ~writeout

        %% plot errors

        figure

        sgtitle(sprintf('Errors, m = %d, n = %d, s = %d, repeats = %d',m,n,sp,num_repeats));

        subplot(1,3,1)

        hold on

        choose_logy = true;

        plot_array = plot_minmax_median_quantiles('-',min_res,max_res,median_res,quant25_res,quant75_res,choose_logy,method_array,iter_save,maxiter,minmaxcolor_dict,quantcolor_dict,linecolor_dict,displayname_dict);

        leg = legend(plot_array, 'location', 'southwest');

        xlabel('$k$','Interpreter','latex')

        ylabel('$$\|Ax-\hat b\|/\|\hat b\|$$','Interpreter','latex')

        hold off





        subplot(1,3,2)

        hold on

        choose_logy = true;

        plot_minmax_median_quantiles('-',min_grad_zfunctional,max_grad_zfunctional,median_grad_zfunctional,...
                                                        quant25_grad_zfunctional,quant75_grad_zfunctional,choose_logy, ...
                                                        method_array,iter_save,maxiter,minmaxcolor_dict,quantcolor_dict,...
                                                        linecolor_dict,displayname_dict);

        xlabel('$k$','Interpreter','latex')
        ylabel('$$\|\nabla_x g^*(\hat b-Ax)\|/\|\hat b\|$$','Interpreter','latex')

        % remove legend
        leg = legend('figure()');
        set(leg,'visible','off')

        hold off









        %% plot results for 'residuals' with P_{R(A)}(b) instead of b
        %{
        subplot(2,3,2)

        hold on

        choose_logy = true;

        medianplot_array = plot_minmax_median_quantiles('-',min_res_proj,max_res_proj,median_res_proj,quant25_res_proj,quant75_res_proj, choose_logy, method_array,iter_save,maxiter,minmaxcolor_dict,quantcolor_dict,linecolor_dict,displayname_dict);

        hold off

        title('Proj. Residuals')
        %}





        %% Plot results for least squares residuals
    %{

        subplot(2,2,3)

        hold on 

        choose_logy = true;

        medianplot_array = plot_minmax_median_quantiles('-',min_lsres,max_lsres,median_lsres,quant25_lsres,quant75_lsres,choose_logy,method_array,iter_save,maxiter,minmaxcolor_dict,quantcolor_dict,linecolor_dict,displayname_dict);

        %legend(medianplot_array, 'location', 'southwest');

        hold off

        title('Gradient of least squares function')
    %}





        %% Plot results for distance to the sparse solution (sol of reg BP problem)


        subplot(1,3,3)

        hold on 

        choose_logy = true;

        plot_minmax_median_quantiles('-',min_err_to_sparse,max_err_to_sparse,median_err_to_sparse,...
                                                        quant25_err_to_sparse,quant75_err_to_sparse, choose_logy,...
                                                        method_array,iter_save,maxiter,minmaxcolor_dict,...
                                                        quantcolor_dict,linecolor_dict,displayname_dict);

        %legend(medianplot_array, 'location', 'northeast');

        % remove legend
        leg = legend('figure()');
        set(leg,'visible','off')

        xlabel('$k$','Interpreter','latex')
        ylabel('$$\|x-\hat x\|/\|\hat x\|$$','Interpreter','latex')

        hold off                                                








        %% Plot results for distance to the Moore-Penrose inverse solution 
        %{
        subplot(1,4,4)

        hold on 

        choose_logy = true;

        plot_minmax_median_quantiles('-',min_err_to_moorepi,max_err_to_moorepi,median_err_to_moorepi,...
                                                        quant25_err_to_moorepi,quant75_err_to_moorepi,choose_logy,...);
                                                        method_array,iter_save,maxiter,minmaxcolor_dict,...
                                                        quantcolor_dict,linecolor_dict,displayname_dict);
        hold off                                                

        title('$$\|x-x^\dagger\|/\|x^\dagger\|$$','Interpreter','latex')
        %}






        %% Plot quantities for stopping criterion

        method_counters_with_stopping = [find(strcmp(method_array,{'grek_1'})), find(strcmp(method_array,{'grek_2'}))...
                                         find(strcmp(method_array,{'grek_3'})), find(strcmp(method_array,{'egrek'}))]; % plot only for these methods
        grek_method_array = {method_array{strcmp(method_array,{'grek_1'})},method_array{strcmp(method_array,{'grek_2'})},...
                             method_array{strcmp(method_array,{'grek_3'})}, method_array{strcmp(method_array,{'egrek'})}};




        figure

        sgtitle(sprintf('Quantities for stopping criterion, solid line: deterministic one, \n dashed line: Stochastic substitute computed in the algorithm'))



        %  plot stochastic substitute for ||A*x_k + b - zstar_k||  together with ||A*x_k + b - zstar_k||
        subplot(1,2,1)
        hold on
        choose_logy = true;

        plot_array = plot_minmax_median_quantiles('-',...
            min_resAbz_mean_list(:,method_counters_with_stopping),...
            max_resAbz_mean_list(:,method_counters_with_stopping),...
            median_resAbz_mean_list(:,method_counters_with_stopping),...
            quant25_resAbz_mean_list(:,method_counters_with_stopping),...
            quant75_resAbz_mean_list(:,method_counters_with_stopping),...
            choose_logy,grek_method_array,iter_save,maxiter,minmaxcolor_dict,quantcolor_dict,linecolor_dict,displayname_dict);

        plot_minmax_median_quantiles(':',...
            min_resAbz_list(:,method_counters_with_stopping),...
            max_resAbz_list(:,method_counters_with_stopping),...
            median_resAbz_list(:,method_counters_with_stopping),...
            quant25_resAbz_list(:,method_counters_with_stopping),...
            quant75_resAbz_list(:,method_counters_with_stopping),...
            choose_logy,grek_method_array,iter_save,maxiter,minmaxcolor_dict,quantcolor_dict,linecolor_dict,displayname_dict);

        if length(method_counters_with_stopping) > 1
           legend(plot_array, 'location', 'northwest');
        end

        xlabel('$$k$$','Interpreter','latex')
        ylabel('$$\|Ax+b-z^*||_2/\|\hat b\|$$','Interpreter','latex')

        hold off



        % plot stochastic substitute for ||ATz_k|| together with ||ATz_k||
        subplot(1,2,2) 
        hold on
        choose_logy = true;

        plot_array = plot_minmax_median_quantiles('-',...
            min_resATz_list(:,method_counters_with_stopping),...
            max_resATz_list(:,method_counters_with_stopping),...
            median_resATz_list(:,method_counters_with_stopping),...
            quant25_resATz_list(:,method_counters_with_stopping),...
            quant75_resATz_list(:,method_counters_with_stopping),...
            choose_logy,grek_method_array,iter_save,maxiter,minmaxcolor_dict,quantcolor_dict,linecolor_dict,displayname_dict);

        plot_minmax_median_quantiles(':',...
            min_resATz_mean_list(:,method_counters_with_stopping),...
            max_resATz_mean_list(:,method_counters_with_stopping),...
            median_resATz_mean_list(:,method_counters_with_stopping),...
            quant25_resATz_mean_list(:,method_counters_with_stopping),...
            quant75_resATz_mean_list(:,method_counters_with_stopping),...
            choose_logy,grek_method_array,iter_save,maxiter,minmaxcolor_dict,quantcolor_dict,linecolor_dict,displayname_dict);

        if length(method_counters_with_stopping) > 1
           legend(plot_array, 'location', 'northwest');
        end

        xlabel('$$k$$','Interpreter','latex')
        ylabel('$$\|A^Tz\|/\|\hat b\|$$','Interpreter','latex')

        hold off

        %% Plot quantities for stopping criterion for only last iterate

        figure

        % plot stochastic substitute for ||Ax-b+z^star||_2^2/norm(b)
        subplot(1,2,1)
        hold on
        choose_logy = true;
        plot_stopped_quantity(resAbz_mean_list,iterstop_list,choose_logy,method_array,method_counter,iter_save,maxiter,linecolor_dict,displayname_dict);
        title('$$Stoch. subst. for \|Ax+b-z^*||_2/\|\hat b\|$$ for last iterate','Interpreter','latex')

        % plot stochastic substitute for ||A^Tz_k||_2^2/norm(b)
        subplot(1,2,2)
        hold on
        choose_logy = true;
        plot_stopped_quantity(resATz_mean_list,iterstop_list,choose_logy,method_array,method_counter,iter_save,maxiter,linecolor_dict,displayname_dict);
        title('$$Stoch. subst. for \|A^Tz||_2/\||\hat b\|$$ for last iterate','Interpreter','latex')








        %% Plot error quantities after stopping for the last iterate x_k for the last instance
          % (only for SREK and ESREK method, each method in a single plot)



        % only plot stopping iteration after

        figure

        sgtitle('Errors for one random instance, circle indicates value at stopped iteration')

        % plot stopped res
        subplot(1,3,1)
        hold on
        choose_logy = true;
        plot_array = plot_stopped_quantity(res,iterstop_list,choose_logy,method_array,method_counter,iter_save,maxiter,linecolor_dict,displayname_dict);
        xlabel('$k$','Interpreter','latex')
        ylabel('$$\|Ax-\hat b\|/\|\hat b\|$$','Interpreter','latex')
        legend(plot_array, 'location', 'northwest');

        % plot stopped grad z functional
        subplot(1,3,2)
        hold on
        choose_logy = true;
        plot_stopped_quantity(grad_zfunctional,iterstop_list,choose_logy,method_array,method_counter,iter_save,maxiter,linecolor_dict,displayname_dict);
        xlabel('$k$','Interpreter','latex')
        ylabel('$$\|A^T \nabla g^*(\hat b-Ax)\|/\|\hat b\|$$','Interpreter','latex')

        % plot stopped reconstruction error 
        subplot(1,3,3)
        hold on
        choose_logy = true;
        plot_stopped_quantity(err_to_sparse,iterstop_list,choose_logy,method_array,method_counter,iter_save,maxiter,linecolor_dict,displayname_dict);
        xlabel('$k$','Interpreter','latex')
        ylabel('$$\|x-\hat x\|/\|\hat x\|$$','Interpreter','latex')





        %% Plot entries of b_exact vs. and b xhat vs. the last iterate x_k 

        % b_exact vs. b

        if real_setting
            figure
            plot_array = plot(1:m,b_exact,1:m,b,':'); title('b_{hat}, b');
        else
            figure
            plot_array = plot(1:m,abs(b_exact),1:m,abs(b),':'); title('|b_{hat}|, |b|');
        end

        legend(plot_array, {'b','noisy b'}, 'location', 'northwest')







        for method_counter = 1:length(method_array)    

            % plot last iterate
            figure
            x = x_last_test_instance(:,method_counter);
            if ~real_setting  % plot absolute values
                xhat = abs(xhat);
                x = abs(x);
            end
            stem(xhat,'color','blue')
            hold on
            stem(x,'color','red')
            title(['After all iterations, ' method_array{method_counter} ' method: x_{hat} (blue), nnz_{hat}=', num2str(sp),', x (red), nnz=', num2str(length(find(abs(x)>1e-5))),' , M-Mult=',num2str(round(iter/max(m,n)))]);
            leg = legend({'x_{hat}','x'}, 'location', 'northwest');

            if any( strcmp( method_array{method_counter}, {'grek_1','grek_2','egrek'} )  )

                % plot iterate after stopping
                figure
                x = xstop_list(:,num_repeats,method_counter);
                if ~real_setting  % plot absolute values
                    xhat = abs(xhat);
                    x = abs(x);
                end        
                stem(xhat,'color','blue'); 
                hold on
                stem(x,'color','red')
                legend({'x_{hat}','x'}, 'location', 'northwest')
                title('After stopping');

            end

        end

    end
    
    
    
    
    
    
    %% generate & save subplot for paper (if writeout) 

    if writeout
        
        figure
  
        % ||Ax-\hat b||
        sgtitle(sprintf('Errors, m = %d, n = %d, s = %d, repeats = %d',m,n,sp,num_repeats));
        subplot(1,3,1)
        choose_logy = true;
        plot_array = plot_minmax_median_quantiles('-',min_res,max_res,median_res,quant25_res,quant75_res,choose_logy,method_array,iter_save,maxiter,minmaxcolor_dict,quantcolor_dict,linecolor_dict,displayname_dict);
        leg = legend(plot_array, 'location', 'southwest');
        if writeout
            set(leg,'visible','off')  % remove legend
        end

        xlabel('$k$','Interpreter','latex')
        ylabel('$$\|Ax-\hat b\|/\|\hat b\|$$','Interpreter','latex')
        


        % ||A^T\nabla g^*(\hat b-Ax)||
        subplot(1,3,2)
        choose_logy = true;
        plot_minmax_median_quantiles('-',min_grad_zfunctional,max_grad_zfunctional,median_grad_zfunctional,...
                                                        quant25_grad_zfunctional,quant75_grad_zfunctional,choose_logy, ...
                                                        method_array,iter_save,maxiter,minmaxcolor_dict,quantcolor_dict,...
                                                        linecolor_dict,displayname_dict);
        xlabel('$k$','Interpreter','latex')
        ylabel('$$\|A^T \nabla g^*(\hat b-Ax)\|/\|\hat b\|$$','Interpreter','latex')

        % remove legend
        leg = legend('figure()');
        set(leg,'visible','off')



        % || x-\hat x||   
        subplot(1,3,3)
        choose_logy = true;
        plot_minmax_median_quantiles('-',min_err_to_sparse,max_err_to_sparse,median_err_to_sparse,...
                                                        quant25_err_to_sparse,quant75_err_to_sparse, choose_logy,...
                                                        method_array,iter_save,maxiter,minmaxcolor_dict,...
                                                        quantcolor_dict,linecolor_dict,displayname_dict);

        %legend(medianplot_array, 'location', 'northeast');

        % remove legend
        leg = legend('figure()');
        set(leg,'visible','off')  
        xlabel('$k$','Interpreter','latex')
        ylabel('$$\|x-\hat x\|/\|\hat x\|$$','Interpreter','latex')   
        
        
        % save first row of figures with matlab2tikz
        matlab2tikz('width','\figurewidth',...
        'extraaxisoptions','legend style={font=\scriptsize},', ...
        [dir_to_figures '/' fig_folder_name '/err_over_iter.tex']);



        % stem plots of last iterates
        
        figure
        
        % b_exact vs. b       

        subplot(1,5,1)

        if real_setting 
            plot(1:m,b_exact,1:m,b,':'); 
        else
            plot(1:m,abs(b_exact),1:m,abs(b),':'); 
        end
        
        axis square

        % remove legend and xlabels
        leg = legend('figure()');
        set(leg,'visible','off')
        set(gca,'XTick',[])
        
         
        plot_nr = 2;

        for method_counter = 1:length(method_array) 

            subplot(1,5,plot_nr)

            % only plot last iterate
            x = x_last_test_instance(:,method_counter);
            
            % prepare x and xhat for getting plot
            tol_plot = 1e-5;
            %xhat = xhat(abs(x) > tol_plot);
            %x = x(abs(x) > tol_plot);
            if ~real_setting  % plot absolute values
                xhat = abs(xhat);
                x = abs(x);
            end
            
            stem(xhat,'color','blue')
            hold on
            stem(x,'color','red')
            axis square
            leg = legend({'x_{hat}','x'}, 'location', 'northwest');
            set(leg,'visible','off')   
            set(gca,'XTick',[])

            plot_nr = plot_nr + 1;

        end




        % save figure with matlab2tikz

        matlab2tikz('width','\figurewidth',...
        'extraaxisoptions','legend style={font=\scriptsize},', ...
        [dir_to_figures '/' fig_folder_name '/x_components.tex']);

    
    
    end
    
        





    %% Display number of nonzero entries of last iterate
    fprintf('\n')
    for i = 1:length(method_array)
        disp(['Entries with absolute value > ' num2str(tol_zero) ', ' method_array{i} ' method (min, median,max): ' num2str(min_nonzero_entries(end,i)) ', ' num2str(median_nonzero_entries(end,i)) ', ' num2str(max_nonzero_entries(end,i)) ])
    end



    %% Display computation time

    min_time = min(time,[],1);
    max_time = max(time,[],1);
    median_time = median(time,1);
    for i = 1:length(method_array)
      disp(['Computation time for ' num2str(maxiter) ' iterations, ' method_array{i} ' method(min,median,max): ' num2str(min_time(i)) ', ' num2str(median_time(i)) ', ' num2str(max_time(i))])
    end


    % ---- until here 




    % output data

    % if more than one repeat, save intersting quantities in a data struct
    
    % first convert cell array to struct (otherwise, problems)
    headings = cell(1,num_methods);
    for i = 1:num_methods
        headings{i} = num2str(i);
    end
    methods = cell2struct(headings,method_array,2);
    
    data = struct('methods', methods, 'iterstop_list', iterstop_list, 'xstop_list', xstop_list, ...
                  'max_res', max_res, 'min_res', min_res, 'median_res', median_res, 'quant25_res', quant25_res, 'quant75_res', quant75_res,...
                  'max_lsres', max_lsres, 'min_lsres', min_lsres, 'median_lsres', median_lsres, 'quant25_lsres', quant25_lsres, 'quant75_lsres', quant75_lsres,...
                  'max_err_to_sparse', max_err_to_sparse, 'min_err_to_sparse', min_err_to_sparse, 'median_err_to_sparse', median_err_to_sparse, 'quant25_err_to_sparse', quant25_err_to_sparse, 'quant75_err_to_sparse', quant75_err_to_sparse,...
                  'err_to_sparse_stopped', err_to_sparse_stopped, 'err_to_sparse_end_div_stopped', err_to_sparse_end_div_stopped, ...
                  'max_err_to_moorepi', max_err_to_moorepi, 'min_err_to_moorepi', min_err_to_moorepi, 'median_err_to_moorepi', median_err_to_moorepi, 'quant25_err_to_moorepi', quant25_err_to_moorepi, 'quant75_err_to_moorepi', quant75_err_to_moorepi,...
                  'max_nonzero_entries', max_nonzero_entries, 'min_nonzero_entries', min_nonzero_entries, 'median_nonzero_entries', median_nonzero_entries, 'quant25_nonzero_entries', quant25_nonzero_entries, 'quant75_nonzero_entries', quant75_nonzero_entries,...
                  'min_rel_cond_A', min_rel_cond_A, 'max_rel_cond_A', max_rel_cond_A);



    % if only one repeat
    %data = struct('method_array', method_array, 'xstop_list', xstop_list, 'zstarstop_list', zstarstop_list, 'res', res, 'lsres', lsres, ...
    %              'err_to_sparse', err_to_sparse, 'err_to_moorepi', err_to_moorepi, 'nonzero_entries', nonzero_entries);                           


    % helper functions

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


