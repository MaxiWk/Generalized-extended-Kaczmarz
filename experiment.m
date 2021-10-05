function data = experiment(n,m,sp,real_setting,lambda_value,T,L_gstar,maxiter,num_repeats,iter_save,rowsamp,colsamp,seed_indices,writeout,disp_instance,savestep,stopcrit_sample_pars,method_array,experiment_description)
%function experimentA_kaczmarz_vs_sparse_kaczmarz(n,m,s,maxiter,num_repeats,rowsamp)
% Perform an experiment to explore usual Kaczmarz vs. sparse Kaczmarz (both
% randomized). It uses an n by m Gaussian random matrix and a sparse right 
% hand size with Gaussian entries. The parameters are:
%  n: number of columns
%  m: number of rows
%  sp: number of nozeros in solution 
%  T: gradient g*. We use
%  g(x) = 1/2 ||x||_2^2 + gamma ||x||_1, 
%  g^*(x) = r_epsilon(x) Huber function, 
%  g^*(x) = r_epsilon(x) + tau/2 ||x||_2^2
%  g^*(x) = mu* r_epsilon(x) + 1/2 ||x||_2^2
%  lambda_value: either a number which is then used as lambda, or 'adapt'
%     (in this case an adaptive value dependent on xhat is used)
%  maxiter: maximum number of Kaczmarz-steps (projections, not sweeps)
%  num_repeats: number of repetitions (samples)
% iter_save: each such number of iterations, a data point is added in the error plot
%  rowsamp: Method to sample the rows. Can be
%    'rownorms_squared': probablity proportional to the squares of row-norms
%    ' uniform': Uniform sampling of rows 
%    ' random_probabilities': Sample a vector of probabilities uniformly at
%        random
%  colsamp: Analogously for column sampling in case of extended Kaczmarz
%  writeout: boolean (if true write data, if false only show plots)
%  savestep: Only each "savesteps"th data point will be stored

% methods: Cell array with methods given as strings: 
% 'rk': Kaczmarz, 'srk': Sparse Randomized Kaczmarz, 'esrk': Exact step Sparse Kaczmarz, 
% 'srek': Sparse Randomized Extended Kaczmarz, 'esrek': Exact step Sparse
% Randomized Extended Kaczmarz

% Some standard parameters for an experiment:
% n = 200;        % Number of columns
% m = round(5*n); % Number of rows
% sp = round(n/8); % Number of nonzeros in solution
% maxiter = 40*m; % Number of iterations
% rowsamp = 'rownorms squared';
% num_repeats = 60;% Number of repeats over random instances
% writeout = false;





disp(experiment_description)

fprintf('Experiment A, m=%d, n=%d, sp=%d, maxiter=%d, num_repeats=%d,\n',m,n,sp,maxiter,num_repeats)




% Initialize storage arrays for results
num_methods = length(method_array); 
[res, lsres, res_proj, lsres_proj, grad_zfunctional, err_to_sparse, err_to_moorepi, nonzero_entries]...
                     = deal( zeros(floor(maxiter/iter_save), num_repeats, num_methods) );
x_last_test_instance = zeros(n,num_methods);
[length_resA_sampled, length_resAT_sampled] = deal(stopcrit_sample_pars.length_resA_sampled,...
                                                             stopcrit_sample_pars.length_resAT_sampled);
resA_sampled = zeros(1, stopcrit_sample_pars.length_resA_sampled);
resAT_sampled = zeros(1, stopcrit_sample_pars.length_resAT_sampled);
[resA_mean, resAT_mean] = deal( zeros(1,maxiter) );

iterstop_list = zeros(num_repeats,num_methods);
xstop_list = zeros(n, num_repeats, num_methods);

time = zeros(num_repeats, num_methods);







% Loop over number of repeats 
fprintf('experiment ')
for repeats = 1:num_repeats 

    fprintf('# %d, ',repeats)
    
    rand('state',repeats);
    randn('state',repeats);

    problem_data = set_up_instance(m,n,sp,real_setting,experiment_description);
    [A,b,b_exact,xhat] = deal( problem_data.A, problem_data.b, problem_data.b_exact, problem_data.xhat );
    [tol_resA, tol_resAT] = deal( problem_data.tol_resA, problem_data.tol_resAT );
    %rankA = rank(A);
    
    if disp_instance 
      disp('')
      disp('xhat (ground truth): ')
      disp(xhat)
      disp('b (ground truth)')
      disp(b_exact)
    end
  
    xmoorepi = pinv(A)*b;
    proj_b = A*xmoorepi;
    
    
    
    
    if lambda_value=='adapt'
        lambda = norm(xhat,1)/(sqrt(0.3*m/nnz(xhat))-1);     % parameter for sparse kaczmarz
    else 
        lambda = lambda_value;
    end
    
    S = @(z) max(abs(z)-lambda,0).*sign(z); % Soft thresholding    

    
    norma = sum(A.^2,2);      % squared norms of rows of A
    norma_col = sum(A.^2,1);  % squared norms of cols of A
    norm_ATb = norm(A'*b);    % for relative least squares residual
    
    
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
   
    
    
    
    for method_counter = 1:length(method_array)
      
        rand('state', seed_indices);
        randn('state', seed_indices);
      
        method = method_array{method_counter};
        iter_save_counter = 1;

        % Initialize method
        x = zeros(n,1);  
        xdual = zeros(n,1);  
        zdual = b;
        z = T(zdual);
        
        % precompute expressions for deterministic algorithms (only for small dimensions)
        
        if any(strcmp(method_array,'det_sek'))
          AAT = A*A'; 
        end
        
        if any(strcmp(method_array,'lin_breg'))
          ATA = A'*A;
          ATb = A'*b;
        end
        
        lin_breg_stepsize = 1/norm(A)^2;
        
        tic;
        
        stopped = false;

          for iter = 1:maxiter

              %if mod(iter, round(maxiter/10)) == 0
              %  disp([num2str(iter) ' iterations'])
              %  disp(['Least-squares residual: ' num2str(norm(A'*(A*x-b)))])
              %endif

              r = samp_row(iter);
              a = A(r,:)';

              s = samp_col(iter); % execute this for all algorithms in order to get the same r indices

              switch method

                % Kaczmarz 
                case 'rk'
                  x = x - (a'*x-b(r))/norma(r)*a;

                % Sparse Kaczmarz 
              case 'srk'
                  xdual = xdual -(a'*x-b(r))/norma(r)*a;
                  x = S(xdual);

                % Sparse Kaczmarz with exact stept
                case 'esrk'
                  [x,xdual] = linesearch_shrinkage(x,xdual,a,b(r),lambda);

                % Extended Kaczmarz  
                case 'rek'
                  c = A(:,s);
                  zdual = zdual - (c'*zdual)/(norma_col(s))*c;
                  x = x -(a'*x-b(r)+zdual(r))/norma(r)*a;

                % Sparse Extended Kaczmarz 
                case 'srek'
                  c = A(:,s);
                  zdual = zdual - (c'*z)/(L_gstar*norma_col(s))*c;
                  z = T(zdual);
                  xdual = xdual -(a'*x-b(r)+zdual(r))/norma(r)*a;
                  x = S(xdual);  

                % Sparse Extended Kaczmarz with exact step
                case 'esrek'
                  c = A(:,s);
                  zdual = zdual - (c'*z)/(L_gstar*norma_col(s))*c;
                  z = T(zdual);
                  [x,xdual] = linesearch_shrinkage(x,xdual,a,b(r)-zdual(r),lambda);

                case 'det_sek'
                  zdual = zdual - lin_breg_stepsize*AAT*zdual;
                  z = T(zdual);
                  xdual = xdual - lin_breg_stepsize*A'*(A*x-b+zdual);
                  x = S(xdual);

              case 'lin_breg'
                  xdual = xdual - lin_breg_stepsize*(ATA*x-ATb);
                  x = S(xdual);

              end 


              % check stopping criterion

              resA_sampled(iter) = a'*x-b(r)+zdual(r);
              first_iter_with_row_samples = max(1, iter-length_resA_sampled+1);
              actual_length_resA_sampled = iter - first_iter_with_row_samples + 1;
              resA_mean(iter) = norm(resA_sampled(first_iter_with_row_samples:iter)) *sqrt(length_resA_sampled/max(1,actual_length_resA_sampled));

              resAT_sampled(iter) = c'*z;
              first_iter_with_col_samples = max(1,iter-length_resAT_sampled+1);
              actual_length_resAT_sampled = iter - first_iter_with_col_samples + 1;
              resAT_mean(iter) = norm(resAT_sampled(first_iter_with_col_samples:iter)) * sqrt(length_resAT_sampled/max(1,actual_length_resAT_sampled));

              if( ~stopped ...
                  && iter > stopcrit_sample_pars.min_possible_iter_for_stopping ...
                  && resA_mean(iter) < tol_resA ...
                  && resAT_mean(iter) < tol_resAT )

                  iterstop_list(repeats,method_counter) = iter;
                  xstop_list(:,repeats,method_counter) = x;
                  stopped = true;                  
                  
                  %break; 

              end



              % save errors 
              if mod(iter, iter_save) == 0  
                %disp(num2str(norm(x-xmoorepi)));  %%%%%%%
                res(iter_save_counter, repeats, method_counter) = norm(A*x-b)/norm(b);
                adaptive_res(iter_save_counter, repeats, method_counter) = norm(A*x-b+zdual)/norm(b);
                %% new idea for measuring lsres: Different for sparse g! Take y^hat = A*x^hat instead of b..
                %lsres(iter_save_counter, repeats, method_counter) = norm(A'*(A*x-A*xhat))/norm_ATb;
                %res_proj(iter_save_counter, repeats, method_counter) = norm(A*x - proj_b)/norm(b);
                lsres(iter_save_counter, repeats, method_counter) = norm(A'*(A*x-b))/norm_ATb;
                %lsres_proj(iter_save_counter, repeats, method_counter) = norm(A'*(A*x - proj_b))/norm_ATb; 
                grad_zfunctional(iter_save_counter, repeats, method_counter) = norm(A'*T(b-A*x))/norm_ATb;
                err_to_moorepi(iter_save_counter, repeats, method_counter) = norm(x-xmoorepi)/norm(xmoorepi);
                err_to_sparse(iter_save_counter, repeats, method_counter) = norm(x-xhat)/norm(xhat);
                tol_zero = 1e-5;   % we count entries with abs value > tol_zero
                nonzero_entries(iter_save_counter, repeats, method_counter) = nnz(abs(x) > tol_zero);        
                iter_save_counter = iter_save_counter + 1;
              end

          end % end for loop over iteration

        time(repeats, method_counter) = toc;

        if ~stopped
            iterstop_list(repeats,method_counter) = maxiter;
            xstop_list(:,repeats,method_counter) = x;
        end
        
        if last_repeat 
            x_last_test_instance(:,method_counter) = x;
        end
    
      end  % end for loop over method

    end  % end for loop over repeats

    disp('')  % new line






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
    linecolor_dict('rk') = 'k';
    linecolor_dict('srk') = 'b';
    linecolor_dict('esrk') = 'y';
    linecolor_dict('rek') = 'c';
    linecolor_dict('srek') = 'r';
    linecolor_dict('esrek') = 'g';
    %linecolor_dict('det_sek') = 'c';
    linecolor_dict('lin_breg') = 'k';

    minmaxcolor_dict = containers.Map();
    minmaxcolor_dict('rk') = lightgray;
    minmaxcolor_dict('srk') = lightblue;
    minmaxcolor_dict('esrk') = lightyellow;
    minmaxcolor_dict('rek') = lightcyan;
    minmaxcolor_dict('srek') = lightred;
    minmaxcolor_dict('esrek') = lightgreen;
    %minmaxcolor_dict('det_sek') = lightcyan;
    minmaxcolor_dict('lin_breg') = lightgray; 

    quantcolor_dict = containers.Map();
    quantcolor_dict('rk') = mediumgray;
    quantcolor_dict('srk') = mediumblue;
    quantcolor_dict('esrk') = mediumyellow;
    quantcolor_dict('rek') = mediumcyan;
    quantcolor_dict('srek') = mediumred; 
    quantcolor_dict('esrek') = mediumgreen;
    %quantcolor_dict('det_sek') = mediumcyan; 
    quantcolor_dict('lin_breg') = mediumgray;

    displayname_dict = containers.Map();
    displayname_dict('rk') = 'Kaczmarz';
    displayname_dict('srk') = 'Sparse Kaczmarz';
    displayname_dict('esrk') = 'Exact step sparse Kaczmarz';
    displayname_dict('rek') = 'Extended Kaczmarz';
    displayname_dict('srek') = 'Sparse extended Kaczmarz'; 
    displayname_dict('esrek') = 'Sparse extended exact step Kaczmarz';
    displayname_dict('det_sek') = 'Deterministic sparse extended Kaczmarz';
    displayname_dict('lin_breg') = 'Linearized Bregman method';


    %% Set up plotable arrays (min, max, quantiles over all repeats)


    % if only one repeat comment out from here ------

    [min_res, max_res, median_res, quant25_res, quant75_res] = compute_minmax_median_quantiles(res);
    [min_adaptive_res, max_adaptive_res, median_adaptive_res, quant25_adaptive_res, quant75_adaptive_res]...
                                                             = compute_minmax_median_quantiles(adaptive_res);
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


    % set zero entries to eps for not getting -inf in log plot
    min_res = set_zero_entries_to_eps(min_res); 
    %min_res_proj = set_zero_entries_to_eps(min_res_proj); 
    min_lsres = set_zero_entries_to_eps(min_lsres); 
    min_grad_zfunctional = set_zero_entries_to_eps(grad_zfunctional);
    %min_lsres_proj = set_zero_entries_to_eps(min_lsres_proj);
    min_adaptive_res = set_zero_entries_to_eps(min_adaptive_res);     
    min_err_to_moorepi = set_zero_entries_to_eps(min_err_to_moorepi);
    min_err_to_sparse = set_zero_entries_to_eps(min_err_to_sparse);





    %% plot results for residuals

    figure

    subplot(2,2,1)

    hold on

    choose_logy = true;

    medianplot_array = plot_minmax_median_quantiles(min_res,max_res,median_res,quant25_res,quant75_res,choose_logy,method_array,iter_save,maxiter,minmaxcolor_dict,quantcolor_dict,linecolor_dict,displayname_dict);

    %legend(medianplot_array, 'location', 'northwest');

    hold off

    title('Residual')




    subplot(2,2,2)

    hold on

    choose_logy = true;

    medianplot_array = plot_minmax_median_quantiles(min_grad_zfunctional,max_grad_zfunctional,median_grad_zfunctional,...
                                                    quant25_grad_zfunctional,quant75_grad_zfunctional,choose_logy, ...
                                                    method_array,iter_save,maxiter,minmaxcolor_dict,quantcolor_dict,...
                                                    linecolor_dict,displayname_dict);


    hold off

    title('Gradient gstar(b-Ax)')




    %% plot results for adaptive residuals

    %{
    subplot(2,2,2)

    hold on

    choose_logy = true;

    medianplot_array = plot_minmax_median_quantiles(min_adaptive_res, max_adaptive_res,...
                    median_adaptive_res, quant25_adaptive_res, quant75_adaptive_res, choose_logy,method_array,iter_save,maxiter,minmaxcolor_dict,quantcolor_dict,linecolor_dict,displayname_dict);

    hold off

    title('||Ax_k-b+z_k*||')
    %}




    %% plot results for 'residuals' with P_{R(A)}(b) instead of b
    %{
    subplot(2,3,2)

    hold on

    choose_logy = true;

    medianplot_array = plot_minmax_median_quantiles(min_res_proj,max_res_proj,median_res_proj,quant25_res_proj,quant75_res_proj, choose_logy, method_array,iter_save,maxiter,minmaxcolor_dict,quantcolor_dict,linecolor_dict,displayname_dict);

    hold off

    title('Proj. Residuals')
    %}





    %% Plot results for least squares residuals


    subplot(2,2,3)

    hold on 

    choose_logy = true;

    medianplot_array = plot_minmax_median_quantiles(min_lsres,max_lsres,median_lsres,quant25_lsres,quant75_lsres,choose_logy,method_array,iter_save,maxiter,minmaxcolor_dict,quantcolor_dict,linecolor_dict,displayname_dict);

    %legend(medianplot_array, 'location', 'southwest');

    hold off

    title('Gradient of least squares function')






    %% Plot results for distance to the sparse solution (sol of reg BP problem)


    subplot(2,2,4)

    hold on 

    choose_logy = true;

    medianplot_array = plot_minmax_median_quantiles(min_err_to_sparse,max_err_to_sparse,median_err_to_sparse,...
                                                    quant25_err_to_sparse,quant75_err_to_sparse, choose_logy,...
                                                    method_array,iter_save,maxiter,minmaxcolor_dict,...
                                                    quantcolor_dict,linecolor_dict,displayname_dict);

    %legend(medianplot_array, 'location', 'northeast');

    hold off                                                

    title('Distance to sparse solution (RegBP)')




    %% Plot results for distance to the sparse solution (sol of reg BP problem)
    %{
    subplot(2,2,4)

    hold on 

    choose_logy = true;

    medianplot_array = plot_minmax_median_quantiles(min_err_to_moorepi,max_err_to_moorepi,median_err_to_moorepi,...
                                                    quant25_err_to_moorepi,quant75_err_to_moorepi,choose_logy,...);
                                                    method_array,iter_save,maxiter,minmaxcolor_dict,...
                                                    quantcolor_dict,linecolor_dict,displayname_dict);
    hold off                                                

    title([sprintf('Distance to Moore Penrose solution; m = %d, n = %d, s = %d, repeats = %d, sampling: ',m,n,sp,num_repeats) rowsamp])
    %}        

    
    
    
    
    %% Plot quantities after stopping for the last iterate x_k
    
    figure 
    
    for i = 1:length(method_array)
        
        % plot stopped res
        subplot(2,2,1)
        plot_array = plot_stopped_quantity(res,iterstop_list,choose_logy,...
                                           method_array,iter_save,maxiter,linecolor_dict,displayname_dict);
        %legend(plot_array, 'location', 'northeast');
        title('Residual')

        % plot stopped grad z functional
        subplot(2,2,2)
        plot_array = plot_stopped_quantity(grad_zfunctional,iterstop_list,choose_logy,...
                                           method_array,iter_save,maxiter,linecolor_dict,displayname_dict);
        %legend(plot_array, 'location', 'northeast');
        title('Gradient gstar(b-Ax)')
        
        % plot stopped lsres
        subplot(2,2,3)
        plot_array = plot_stopped_quantity(lsres,iterstop_list,choose_logy,...
                                           method_array,iter_save,maxiter,linecolor_dict,displayname_dict);
        %legend(plot_array, 'location', 'northeast');
        title('Gradient of least squares function')
        
        % plot stopped reconstruction error 
        subplot(2,2,4)
        plot_array = plot_stopped_quantity(err_to_sparse,iterstop_list,choose_logy,...
                                           method_array,iter_save,maxiter,linecolor_dict,displayname_dict);
        %legend(plot_array, 'location', 'northeast');        
        title('Distance to sparse solution (RegBP)')
        
    end


    %% Plot entries of b_exact vs. and b xhat vs. the last iterate x_k 
    
    % b_exact vs. b
    
    if real_setting
        figure
        plot(1:m,b_exact,1:m,b,':'); title('b_{hat}, b');
        title('Last test instance: hidden b in R(A) and noisy b')
    else
        figure
        plot(1:m,abs(b_exact),1:m,abs(b),':'); title('|b_{hat}|, |b|');
        title('Last test instance: abs. value of b in R(A) and noisy b')
    end
        
    
    for method_counter = 1:length(method_array)           
        
        % plot last iterate
        figure
        x = x_last_test_instance(:,method_counter);
        if ~real_setting  % plot absolute values
            xhat = abs(xhat);
            x = abs(x);
        end
        stem(xhat,'color','blue'); 
        hold on 
        stem(x,'color','red');
        title(['After all iterations: x_{hat} (blue), nnz_{hat}=', num2str(sp),' , x (red), nnz=', num2str(length(find(abs(x)>1e-5))),' , M-Mult=',num2str(round(iter/max(m,n)))]);
        
        % plot iterate after stopping
        figure
        x = xstop_list(:,num_repeats,method_counter);
        if ~real_setting  % plot absolute values
            xhat = abs(xhat);
            x = abs(x);
        end        
        stem(xhat,'color','blue'); 
        hold on
        stem(x,'color','red');
        title(['After stopping: x_{hat} (blue) , nnz_{hat}=', num2str(sp),' , x (red) , nnz=', num2str(length(find(abs(x)>1e-5))),' , M-Mult=',num2str(round(iter/max(m,n)))]);

    end




    %% Display number of nonzero entries of last iterate

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
                  'max_err_to_moorepi', max_err_to_moorepi, 'min_err_to_moorepi', min_err_to_moorepi, 'median_err_to_moorepi', median_err_to_moorepi, 'quant25_err_to_moorepi', quant25_err_to_moorepi, 'quant75_err_to_moorepi', quant75_err_to_moorepi,...
                  'max_nonzero_entries', max_nonzero_entries, 'min_nonzero_entries', min_nonzero_entries, 'median_nonzero_entries', median_nonzero_entries, 'quant25_nonzero_entries', quant25_nonzero_entries, 'quant75_nonzero_entries', quant75_nonzero_entries);



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


