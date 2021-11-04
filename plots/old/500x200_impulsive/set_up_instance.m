function problem_data = set_up_instance(m,n,sp,real_setting,experiment_description)
  
  b_exact = 0;  
  tol_resA = 0;
  tol_resAT = 0;
  
  switch experiment_description
    
    
    
    
     case 'full rank, solvable, no noise'
       
          A = randn(m,n);           % Operator
          rel_cond_A = compute_rel_cond(A,min(m,n));
          xhat = sparserandn(n,sp);  % true solution
          b_exact = A*xhat;               % exact data
          b = b_exact;
          
          tol_resAbz = 1e-7 * norm(b);
          tol_resATz = 1e-7 * norm(b); 

            
     
        
     case 'rank-deficient, only noise in R(A) complement'
       
          % m<=n and least-squares solution shall be not unique (no full rank)          
          
          noiselev_rangeA_ortho = 0.5;
          rank = round(min(m,n)/2);   % must be smaller than m = min(m,n)
          real_setting = true;
          sing_values_data.distribution = 'normal';
          sing_values_data.mu = 0; sing_values_data.stddev = 1;
          
          A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
                 
          xhat = sparserandn(n,sp);  % true solution          
          b_exact = A*xhat;               % exact data
          
          % create offset with noise in the complement of R(A)
          N = null(A');             % null(A') = orthogonal complement of R(A)
          v = randn(size(N,2),1);
          noise_rangeA_ortho = N*v; 
          noise_rangeA_ortho = noiselev_rangeA_ortho* noise_rangeA_ortho/norm(noise_rangeA_ortho);
          b = b_exact + noise_rangeA_ortho;    
          
          tol_resAbz = 1e-4 *norm(b);
          tol_resATz = 1e-4 *norm(b);
          
          
   
          
      case  'rank deficient, noise split into R(A) and R(A) complement'
  
          % m<=n and least-squares solution shall be not unique (no full rank), finally also add noise to b         
          
          noiselev_rangeA_ortho = 0.5;
          noiselev_rangeA = 0.1;
          rank = round(min(m,n)/2);   % must be smaller than m = min(m,n)
          
          real_setting = true;
          sing_values_data.distribution = 'normal';
          sing_values_data.mu = 0; sing_values_data.stddev = 1;
          
          A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
                  
          xhat = sparserandn(n,sp);  % true solution          
          b_exact = A*xhat;               % exact data
          
          % create noisy b 

          N = null(A');             % null(A') = orthogonal complement of R(A)
          v = randn(size(N,2),1);
          rangeA_ortho = N*v;
          rangeA_ortho = noiselev_rangeA_ortho* rangeA_ortho/norm(rangeA_ortho);
          noise_rangeA = A* randn(n,1);
          noise_rangeA = noiselev_rangeA/ norm(noise_rangeA) *noise_rangeA;

          b = b_exact + rangeA_ortho + noise_rangeA;      

          tol_resAbz = 0.1*noiselev_rangeA *norm(b);
          tol_resATz = 0.01;
      
          
          
          
          
    case 'rank deficient, only noise in R(A) complement, well conditioned A'
      
      noiselev_rangeA_ortho = 0.5;
      rank = round(min(m,n)/2);
      
      sing_values_data.distribution = 'uniform';
      sing_values_data.sigma_min = 1; sing_values_data.sigma_max = 2;

      A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
      rel_cond_A = compute_rel_cond(A,rank);
          
      xhat = sparserandn(n,sp);  % true solution
      b_exact = A*xhat;               % exact data
      N = null(A');             % null(A') = orthogonal complement of R(A)
      v = randn(size(N,2),1);
      noise_rangeA_ortho = N*v;
      noise_rangeA_ortho = noiselev_rangeA_ortho* noise_rangeA_ortho/norm(noise_rangeA_ortho);
      b = b_exact + noise_rangeA_ortho;

      tol_resAbz = 1e-14*norm(b);
      tol_resATz = 1e-14*norm(b);
          
 
    
    case 'rank deficient, only noise in R(A) complement, well conditioned A and xhat'  
        
      noiselev_rangeA_ortho = 5;
      sing_values_data.distribution = 'uniform';
      sing_values_data.sigma_min = 1;
      sing_values_data.sigma_max = 2;
      xhat_min = 1;
      xhat_max = 6;
      rank = round(min(m,n)/2);

      A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
      rel_cond_A = compute_rel_cond(A,rank);
               
      % create xhat

      xhat=zeros(n,1);
      rpn = randperm(n); 
      nnz_ind = rpn(1:sp);
      if real_setting
          xhat(nnz_ind)= xhat_min+ (xhat_max-xhat_min)*rand(sp,1);
      else
          xhat(nnz_ind)= xhat_min+ (xhat_max-xhat_min)*(rand(sp,1) +1i*rand(sp,1));
      end
        
      % create noisy b 
      
      b_exact = A*xhat;               % exact data
      N = null(A');             % null(A') = orthogonal complement of R(A)
      v = randn(size(N,2),1);
      rangeA_ortho = N*v;
      rangeA_ortho = noiselev_rangeA_ortho* rangeA_ortho/norm(rangeA_ortho);
      
      b = b_exact + rangeA_ortho;

      tol_resAbz = 1e-12*norm(b);
      tol_resATz = 1e-12*norm(b);
          
          
      
      
      
    case 'rank deficient, noise split into R(A) and R(A) complement, well conditioned A and xhat'  
        
      noiselev_rangeA_ortho = 5;
      noiselev_rangeA = 0.01;
      sing_values_data.distribution = 'uniform';
      sing_values_data.sigma_min = 1;
      sing_values_data.sigma_max = 2;
      xhat_min = 1;
      xhat_max = 6;
      rank = round(min(m,n)/2);  

      A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
      rel_cond_A = compute_rel_cond(A,rank);
                          
      % create xhat
      xhat = random_xhat_uniform_from_sp_min_max(sp, xhat_min, xhat_max);

        
      % create noisy b 
      
      b_exact = A*xhat;               % exact data
      N = null(A');             % null(A') = orthogonal complement of R(A)
      v = randn(size(N,2),1);
      rangeA_ortho = N*v;
      rangeA_ortho = noiselev_rangeA_ortho* rangeA_ortho/norm(rangeA_ortho);
      noise_rangeA = A* randn(n,1);
      noise_rangeA = noiselev_rangeA/ norm(noise_rangeA) *noise_rangeA;
      
      b = b_exact + rangeA_ortho + noise_rangeA;      
          
      tol_resAbz = noiselev_rangeA *norm(b);i
      tol_resATz = 0.1;
      
      
      
      
     case 'dctmatrix, noise split into R(A) and R(A) complement'
         
       noiselev_rangeA_ortho = 5;
       noiselev_rangeA = 0.1;  % auch mal 0.1
       rank = 500;
         
       U=randn(m,rank); V=[eye(rank);dctmtx(rank)];
       A=U*V';
       rel_cond_A = compute_rel_cond(A,rank);
       
       xhat_min = 1;
       xhat_max = 6;
       xhat = random_xhat_uniform_from_sp_min_max(sp, xhat_min, xhat_max);
       
       % create noisy b 
      
       b_exact = A*xhat;               % exact data
       N = null(A');             % null(A') = orthogonal complement of R(A)
       v = randn(size(N,2),1);
       rangeA_ortho = N*v;
       rangeA_ortho = noiselev_rangeA_ortho* rangeA_ortho/norm(rangeA_ortho);
       noise_rangeA = A* randn(n,1);
       noise_rangeA = noiselev_rangeA/ norm(noise_rangeA) *noise_rangeA;
      
       b = b_exact + rangeA_ortho + noise_rangeA;      
          
       tol_resAbz = noiselev_rangeA *norm(b);
       tol_resATz = 0.1;
      
       
       
       
       
      case 'full rank, small gaussian noise'
   
          noiselev = 0.1;
     
          A = randn(m,n);           % Operator
          rel_cond_A = compute_rel_cond(A,min(m,n));
          xhat = sparserandn(n,sp);  % true solution
          b_exact = A*xhat;               % exact data
          
          noise = randn(size(b_exact));
          b = b_exact + noiselev*noise/norm(noise);
          
          tol_resAbz = 1e-7 * norm(b);
          tol_resATz = 0.1;  
          
          
          
    case 'full rank, medium gaussian noise'
      
          noiselev = 0.5;
          
          A = randn(m,n);           % Operator
          rel_cond_A = compute_rel_cond(A,min(m,n));
          xhat = sparserandn(n,sp);  % true solution
          b_exact = A*xhat;               % exact data
          
          noise = randn(size(b_exact));
          b = b_exact + noiselev*noise/norm(noise);

          tol_resAbz = 1e-7 * norm(b);
          tol_resATz = 0.5;  
          
          
          
     case 'full rank, large gaussian noise'
      
          noiselev = 2;
          
          A = randn(m,n);           % Operator
          rel_cond_A = compute_rel_cond(A,min(m,n));
          xhat = sparserandn(n,sp);  % true solution
          b_exact = A*xhat;               % exact data
          
          noise = randn(size(b_exact));
          b = b_exact + noiselev*noise/norm(noise);
          
          tol_resAbz = 1e-7 * norm(b);
          tol_resATz = 0.5;            
          
          
          
    case 'large uniform noise, full rank'
      
          noiselev = 0.5;
          
          A = randn(m,n);           % Operator
          rel_cond_A = compute_rel_cond(A,min(m,n));
          xhat = sparserandn(n,sp);  % true solution
          b_exact = A*xhat;               % exact data
          
          noise = noiselev* (1 - 2*rand(size(b_exact)));
          b = b_exact + noise;

          tol_resAbz = 1e-7 * norm(b);
          tol_resATz = 0.5;  
          
          


    case 'rank deficient, small uniform noise'
      
          noiselev = 0.01;
          rank = round(min(m,n)/2); 
          sing_values_data.distribution = 'normal';
          sing_values_data.mu = 0;
          sing_values_data.stddev = 1;
          A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
          xhat = sparserandn(n,sp);
          b_exact = A*xhat;
          
          noise = noiselev* (1 - 2*rand(size(b_exact)));
          b = b_exact + noise;         

          tol_resAbz = 1e-7 * norm(b);
          tol_resATz = 0.1;   




    case 'rank deficient, small gaussian noise'
      
          noiselev = 0.01;
          rank = round(min(m,n)/2); 
          real_setting = true;
          sing_values_data.distribution = 'normal';
          sing_values_data.mu = 0; 
          sing_values_data.stddev = 1;
          A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
          xhat = sparserandn(n,sp);
          b_exact = A*xhat;
          
          noise = noiselev* randn(size(b_exact));
          b = b_exact + noise;   
         
          tol_resAbz = 1e-7 * norm(b);
          tol_resATz = 0.5; 
          
          
          
        
    case 'rank deficient, medium uniform noise'
    
          noiselev = 0.1;
          rank = round(min(m,n)/2);
          real_setting = true;
          sing_values_data.distribution = 'normal';
          sing_values_data.mu = 0; 
          sing_values_data.stddev = 1;
          A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
          xhat = sparserandn(n,sp);
          b_exact = A*xhat;
          
          noise = noiselev* (1 - 2*rand(size(b_exact)));
          b = b_exact + noise;

          tol_resAbz = 1e-7 * norm(b);
          tol_resATz = 0.5; 

          
          
  case 'impulsive noise, full rank'
  
          num_comp_noise = 10;
          noiselev_impulsive_noise_factor = 0.1;

          A = randn(m,n);                  
          rel_cond_A = compute_rel_cond(A,min(m,n));
          xhat = sparserandn(n,sp);  % true solution          
          b_exact = A*xhat;               % exact data
          
          num_comp_noise = min(m,num_comp_noise);
          perm = randperm(m);
          b = b_exact;
          
          noiselev_impulsive = noiselev_impulsive_noise_factor *norm(b);
          noise = randn(num_comp_noise,1);
          b(perm(1:num_comp_noise)) = b(perm(1:num_comp_noise)) + noiselev_impulsive* noise / norm(noise); 
          
          tol_resAbz = noiselev_impulsive_noise_factor *noiselev_impulsive;
          tol_resATz = 1e-5; 
 
 
          
  case 'impulsive noise, rank deficient'
  
          num_comp_noise = 10;
          noiselev_impulsive_noise_factor = 10;
          
          rank = round(min(m,n)/2);
          real_setting = true;
          sing_values_data.distribution = 'normal';
          sing_values_data.mu = 0;
          sing_values_data.stddev = 1;

          A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
          xhat = sparserandn(n,sp);  % true solution          
          b_exact = A*xhat;               % exact data
          
          num_comp_noise = min(m,num_comp_noise);
          perm = randperm(m);
          b = b_exact;
          
          noiselev_impulsive = noiselev_impulsive_noise_factor *norm(b);
          noise = randn(num_comp_noise,1);
          b(perm(1:num_comp_noise)) = b(perm(1:num_comp_noise)) + noiselev_impulsive* noise / norm(noise); 
          
          tol_resAbz = noiselev_impulsive_noise_factor*noiselev_impulsive*1e-6;
          tol_resATz = noiselev_impulsive_noise_factor*noiselev_impulsive*1e-4; 
          
    
          
  case 'impulsive noise, rank-deficient, A with unif distr sv'
  
          num_comp_noise = 10;
          noiselev_impulsive_noise_factor = 0.1;
          
          rank = round(min(m,n)/2);
          sing_values_data.distribution = 'uniform';
          sing_values_data.sigma_min = 0.001;
          sing_values_data.sigma_max = 1;

          A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
                    
          xhat = sparserandn(n,sp);  % true solution          
          b_exact = A*xhat;               % exact data
          
          num_comp_noise = min(m,num_comp_noise);
          perm = randperm(m);
          b = b_exact;
          
          noiselev_impulsive = noiselev_impulsive_noise_factor *norm(b);
          noise = randn(num_comp_noise,1);
          b(perm(1:num_comp_noise)) = b(perm(1:num_comp_noise)) + noiselev_impulsive* noise / norm(noise); 
          
          tol_resAbz = noiselev_impulsive_noise_factor *noiselev_impulsive;
          tol_resATz = 1e-5; 
          
          
          
          
    case 'impulsive noise, rank-deficient, well conditioned A and xhat'
  
          num_comp_noise = 10;
          noiselev_impulsive_noise_factor = 10;

          rank = round(min(m,n)/2);
          sing_values_data.distribution = 'uniform';
          sing_values_data.sigma_min = 1;
          sing_values_data.sigma_max = 2;

          A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
                    
          xhat = sparserandn(n,sp);  % true solution          
          b_exact = A*xhat;               % exact data
          
          num_comp_noise = min(m,num_comp_noise);
          perm = randperm(m);
          b = b_exact;
          
          noiselev_impulsive = noiselev_impulsive_noise_factor *norm(b);
          noise = randn(num_comp_noise,1);
          b(perm(1:num_comp_noise)) = b(perm(1:num_comp_noise)) + noiselev_impulsive* noise / norm(noise); 
          
          tol_resAbz = 1e-4 *noiselev_impulsive_noise_factor *noiselev_impulsive;
          tol_resATz = 1e-5;         

          
       
   case 'impulsive noise and additional noise, rank-deficient, well conditioned A and xhat'
  
          num_comp_noise = 10;
          noiselev_impulsive_noise_factor = 0.1;
          noiselev_add_noise = 0.01;

          rank = round(min(m,n)/2);
          sing_values_data.distribution = 'uniform';
          sing_values_data.sigma_min = 1;
          sing_values_data.sigma_max = 2;

          A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
                    
          xhat = sparserandn(n,sp);  % true solution          
          b_exact = A*xhat;          % exact data
          
          num_comp_noise = min(m,num_comp_noise);
          perm = randperm(m);
          b = b_exact;
          
          noiselev_impulsive = noiselev_impulsive_noise_factor*norm(b);
          noise = randn(num_comp_noise,1);
          b(perm(1:num_comp_noise)) = b(perm(1:num_comp_noise)) + noiselev_impulsive* noise / norm(noise); 
          
          tol_resAbz = noiselev_impulsive_noise_factor *noiselev_impulsive;
          tol_resATz = 1e-5;              
 
          
    
     case 'impulsive noise Franks example' % m=800, n=600, rank=600
         
          num_comp_noise = ceil(m/40);
          noiselev_impulsive_noise_factor = 40;
          noiselev_gaussian_noise_factor = 0.01;
          
          rank = min(m,n);
          sing_values_data.distribution = 'uniform';
          sing_values_data.sigma_min = 1;
          sing_values_data.sigma_max = 2;
          A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
          
          xhat = sparserandn(n,sp);  % true solution          
          b_exact = A*xhat;          % exact data         
          
          perm = randperm(m);
          noise_indices = perm(1:num_comp_noise);

          if real_setting
              b=b_exact; 
              b(noise_indices)=b(noise_indices)...
                  + noiselev_impulsive_noise_factor *max(abs(b))*(rand(length(noise_indices),1)-0.5);
          else
              b=b_hat; 
              b(noise_indices)=b(noise_indices)...
                         + noiselev_impulsive_noise_factor *max(abs(b)) ...
                         * ((rand(noise_indices,1)-0.5)+1i*(rand(length(noise_indices),1)-0.5));
          end
          
          
          % add gaussian noise
          b_add = randn(size(b));
          b = b + noiselev_gaussian_noise_factor * norm(b) * b_add/norm(b_add);
         
          tol_resAbz = noiselev_gaussian_noise_factor * norm(b);
          tol_resATz = 0.1;
          
          
          
     otherwise 
          disp('No experiment chosen!')
 
          
         
  end
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%

  function A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data)
      
      if real_setting
          U=orth(randn(m,rank)); V=orth(randn(n,rank));
      else
          U=orth(randn(m,rank)+1i*randn(m,rank)); V=orth(randn(n,rank)+1i*randn(n,rank));
      end
      
      if strcmp(sing_values_data.distribution,'uniform')
          sing_values = sing_values_data.sigma_min + ...
               (sing_values_data.sigma_max - sing_values_data.sigma_min) *rand(1,rank);
      elseif strcmp(sing_values_data.distribution,'normal')
          sing_values = sing_values_data.mu + sing_values_data.stddev *randn(1,rank);
      else
          error('Failed to create system matrix: No matching string!')
      end
      
      D = diag(sing_values);  
      A=U*D*V';
      
  end
      
      
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

  function xhat = random_xhat_uniform_from_sp_min_max(sp, xhat_min, xhat_max)
      xhat=zeros(n,1);
      rpn = randperm(n); 
      nnz_ind = rpn(1:sp);
      if real_setting
          xhat(nnz_ind)= xhat_min+ (xhat_max-xhat_min)*rand(sp,1);
      else
          xhat(nnz_ind)= xhat_min+ (xhat_max-xhat_min)*(rand(sp,1) +1i*rand(sp,1));
      end
      
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function rel_cond_A = compute_rel_cond(A,rank)
      s = svd(A);
      s = s(1:rank);
      if min(s) < eps
          rel_cond_A = cond(A);
      else
          rel_cond_A = max(s)/min(s);
      end
  end






 problem_data.A = A; 
 problem_data.rel_cond_A = rel_cond_A;
 problem_data.xhat = xhat;
 problem_data.b_exact = b_exact;
 problem_data.b = b;
 problem_data.tol_resAbz = tol_resAbz;
 problem_data.tol_resATz = tol_resATz;

end
