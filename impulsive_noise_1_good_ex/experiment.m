function data = experiment(n,m,sp,lambda_value,T,L_gstar,maxiter,num_repeats,iter_save,rowsamp,colsamp,seed_indices,writeout,disp_instance,savestep,method_array,experiment_description)
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




% Initialize storage for results
num_methods = length(method_array); 
[res, lsres, res_proj, lsres_proj, err_to_sparse, err_to_moorepi, nonzero_entries]...
                     = deal( zeros(floor(maxiter/iter_save), num_repeats, num_methods) );
xhat_list = zeros(n, num_repeats, num_methods);
yhat_list = zeros(m, num_repeats, num_methods); 
time = zeros(num_repeats, num_methods);

fprintf('experiment ')


% Loop over number of repeats 
for repeats = 1:num_repeats 

    fprintf('# %d, ',repeats)
    
    rand('state',repeats);
    randn('state',repeats);

    [A,b,b_exact,xhat] = set_up_instance(m,n,sp,experiment_description);
    rankA = rank(A);
    
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
        lambda = norm(xhat,1)/(sqrt(0.3*m/nnz(xhat))-1)     % parameter for sparse kaczmarz
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
              zdual = zdual - (c'*z)/(L*norma_col(s))*c;
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
          

          if mod(iter, iter_save) == 0  
            %disp(num2str(norm(x-xmoorepi)));  %%%%%%%
            res(iter_save_counter, repeats, method_counter) = norm(A*x-b)/norm(b);
            adaptive_res(iter_save_counter, repeats, method_counter) = norm(A*x-b+zdual)/norm(b);
            %% new idea for measuring lsres: Different for sparse g! Take y^hat = A*x^hat instead of b..
            %lsres(iter_save_counter, repeats, method_counter) = norm(A'*(A*x-A*xhat))/norm_ATb;
            %res_proj(iter_save_counter, repeats, method_counter) = norm(A*x - proj_b)/norm(b);
            lsres(iter_save_counter, repeats, method_counter) = norm(A'*(A*x-b))/norm_ATb;
            %lsres_proj(iter_save_counter, repeats, method_counter) = norm(A'*(A*x - proj_b))/norm_ATb; 
            err_to_moorepi(iter_save_counter, repeats, method_counter) = norm(x-xmoorepi)/norm(xmoorepi);
            err_to_sparse(iter_save_counter, repeats, method_counter) = norm(x-xhat)/norm(xhat);
            tol_zero = 1e-5;   % we count entries with abs value > tol_zero
            nonzero_entries(iter_save_counter, repeats, method_counter) = nnz(abs(x) > tol_zero);        
            iter_save_counter = iter_save_counter + 1;
          end
          
      end % end for loop over iteration
      
    time(repeats, method_counter) = toc;
    
    
    xhat_list(:,repeats,method_counter) = x;
    if any( strcmp(method, {'rek','srek','esrek'}) )
      yhat_list(:,repeats,method_counter) = b-zdual;
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

% set zero entries to eps for not getting -inf in log plot
min_res = set_zero_entries_to_eps(min_res); 
%min_res_proj = set_zero_entries_to_eps(min_res_proj); 
min_lsres = set_zero_entries_to_eps(min_lsres); 
%min_lsres_proj = set_zero_entries_to_eps(min_lsres_proj);
min_adaptive_res = set_zero_entries_to_eps(min_adaptive_res);     
min_err_to_moorepi = set_zero_entries_to_eps(min_err_to_moorepi);
min_err_to_sparse = set_zero_entries_to_eps(min_err_to_sparse);



%% plot results for residuals

figure

subplot(1,4,1)

hold on

choose_logy = true;

medianplot_array = plot_minmax_median_quantiles(min_res,max_res,median_res,quant25_res,quant75_res,choose_logy);

hold off

title([sprintf('Residuals; m = %d, n = %d, rank = %d s = %d, repeats = %d, sampling: ',m,n,rankA,sp,num_repeats) rowsamp])




%% plot results for adaptive residuals

subplot(1,4,2)

hold on

choose_logy = true;

medianplot_array = plot_minmax_median_quantiles(min_adaptive_res, max_adaptive_res,...
                median_adaptive_res, quant25_adaptive_res, quant75_adaptive_res, choose_logy);

hold off

title([sprintf('||Ax_k-b+z_k*||; m = %d, n = %d, s = %d, repeats = %d, sampling: ',m,n,sp,num_repeats) rowsamp])





%% plot results for 'residuals' with P_{R(A)}(b) instead of b
%{
subplot(2,3,2)

hold on

choose_logy = true;

medianplot_array = plot_minmax_median_quantiles(min_res_proj,max_res_proj,median_res_proj,quant25_res_proj,quant75_res_proj, choose_logy);

hold off

title([sprintf('Proj. Residuals; m = %d, n = %d, s = %d, repeats = %d, sampling: ',m,n,sp,num_repeats) rowsamp])
%}




      
%% Plot results for least squares residuals

subplot(1,4,3)

hold on 

choose_logy = true;

medianplot_array = plot_minmax_median_quantiles(min_lsres,max_lsres,median_lsres,quant25_lsres,quant75_lsres,choose_logy);

hold off

title([sprintf('Gradient of least squares function; m = %d, n = %d, rank = %d, s = %d, repeats = %d, sampling: ',m,n,rankA,sp,num_repeats) rowsamp])






%% Plot results for distance to the sparse solution (sol of reg BP problem)


subplot(1,4,4)

hold on 

choose_logy = true;

medianplot_array = plot_minmax_median_quantiles(min_err_to_sparse,max_err_to_sparse,median_err_to_sparse,...
                                                quant25_err_to_sparse,quant75_err_to_sparse, choose_logy);
                                                
hold off                                                

title([sprintf('Distance to sparse solution (RegBP); m = %d, n = %d, rank = %d, s = %d, repeats = %d, sampling: ',m,n,rankA,sp,num_repeats) rowsamp])
                            
                                      
                                 
                             
%% Plot results for distance to the sparse solution (sol of reg BP problem)
%{
subplot(2,3,4)

hold on 

choose_logy = true;

medianplot_array = plot_minmax_median_quantiles(min_err_to_moorepi,max_err_to_moorepi,median_err_to_moorepi,...
                                                quant25_err_to_moorepi,quant75_err_to_moorepi,choose_logy);
                                                
hold off                                                

title([sprintf('Distance to Moore Penrose solution; m = %d, n = %d, s = %d, repeats = %d, sampling: ',m,n,sp,num_repeats) rowsamp])
       %}   
        
      







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

% if more than one repeat

data = struct('method_array', method_array, 'xhat_list', xhat_list, 'yhat_list', yhat_list, ...
              'max_res', max_res, 'min_res', min_res, 'median_res', median_res, 'quant25_res', quant25_res, 'quant75_res', quant75_res,...
              'max_lsres', max_lsres, 'min_lsres', min_lsres, 'median_lsres', median_lsres, 'quant25_lsres', quant25_lsres, 'quant75_lsres', quant75_lsres,...
              'max_err_to_sparse', max_err_to_sparse, 'min_err_to_sparse', min_err_to_sparse, 'median_err_to_sparse', median_err_to_sparse, 'quant25_err_to_sparse', quant25_err_to_sparse, 'quant75_err_to_sparse', quant75_err_to_sparse,...
              'max_err_to_moorepi', max_err_to_moorepi, 'min_err_to_moorepi', min_err_to_moorepi, 'median_err_to_moorepi', median_err_to_moorepi, 'quant25_err_to_moorepi', quant25_err_to_moorepi, 'quant75_err_to_moorepi', quant75_err_to_moorepi,...
              'max_nonzero_entries', max_nonzero_entries, 'min_nonzero_entries', min_nonzero_entries, 'median_nonzero_entries', median_nonzero_entries, 'quant25_nonzero_entries', quant25_nonzero_entries, 'quant75_nonzero_entries', quant75_nonzero_entries);



% if only one repeat
%data = struct('method_array', method_array, 'xhat_list', xhat_list, 'yhat_list', yhat_list, 'res', res, 'lsres', lsres, ...
%              'err_to_sparse', err_to_sparse, 'err_to_moorepi', err_to_moorepi, 'nonzero_entries', nonzero_entries);                           


% helper functions

% for logarithmic plot
function arr = set_zero_entries_to_eps(arr)

    arr(arr<eps) = eps;
    
end




function [mins, maxs, medians, quant25s, quant75s] = compute_minmax_median_quantiles(arr)
    for i = 1:size(arr,3)  % loop over all methods
      mins(:,i) = min(arr(:,:,i), [], 2);
      maxs(:,i) = max(arr(:,:,i), [], 2);
      medians(:,i) = median(arr(:,:,i),2);
      quant25s(:,i) = quantile(arr(:,:,i),0.25,2); 
      quant75s(:,i) = quantile(arr(:,:,i),0.75,2);
    end
end




function medianplot_array = plot_minmax_median_quantiles(mins,maxs,medians,quant25s,quant75s,choose_logy)
  medianplot_array = zeros(1, length(method_array));  % for legend
  num_iter_array = 1:iter_save:maxiter;
  for i = 1:num_methods
      minmaxcolor_i = minmaxcolor_dict(method_array{i});
      
      if choose_logy
        h = fill([num_iter_array  fliplr(num_iter_array)], [log10(maxs(:,i)')  fliplr(log10(mins(:,i))')], minmaxcolor_i,'EdgeColor', 'none');
        set(h,'facealpha', .5)
        quantcolor_i = quantcolor_dict(method_array{i});
        h = fill([num_iter_array  fliplr(num_iter_array)], [log10(quant75s(:,i)')  fliplr(log10(quant25s(:,i))')], quantcolor_i,'EdgeColor', 'none');
        set(h,'facealpha', .5)
        medianplot_array(i) = plot( num_iter_array,log10(medians(:,i)),linecolor_dict(method_array{i}),'LineWidth',2,...
                   'DisplayName',displayname_dict(method_array{i}) );
        ylabel('(log scale)')
      else 
        h = fill([num_iter_array  fliplr(num_iter_array)], [maxs(:,i)'  fliplr(mins(:,i)')], minmaxcolor_i,'EdgeColor', 'none');
        set(h,'facealpha', .5)
        quantcolor_i = quantcolor_dict(method_array{i});
        h = fill([num_iter_array  fliplr(num_iter_array)], [quant75s(:,i)'  fliplr(quant25s(:,i)')], quantcolor_i,'EdgeColor', 'none');
        set(h,'facealpha', .5)
        medianplot_array(i) = plot( num_iter_array,medians(:,i),linecolor_dict(method_array{i}),'LineWidth',2,...
                   'DisplayName',displayname_dict(method_array{i}) );          
      end
      % use 10^ notation for iterations    
      %xt = get(gca, 'xtick');
      %set('xticklabel', sprintf('%1.1e|', xt));
      
  end
  legend(medianplot_array, 'location', 'southwest');
end


end
