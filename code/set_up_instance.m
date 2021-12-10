function problem_data = set_up_instance(m,n,sp,real_setting,experiment_description)
  
  b_exact = 0;  
  tol_resA = 0;
  tol_resAT = 0;
  
  switch experiment_description
    
    
      
    
      
     case 'small_example'
         
         A = [1,2;3,4];
         xhat = [5,1]';
         b = A*xhat;
         
         b_exact = b;
         tol_resAbz = 1e-7;
         tol_resATz = 1e-7;         
         rel_cond_A = cond(A);
         
          
          
    
     case 'full rank, consistent, no noise'
       
          rank = min(m,n);
          sing_values_data.distribution = 'normal';
          sing_values_data.mu = 0; sing_values_data.stddev = 1;
          A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
          rel_cond_A = compute_rel_cond(A,min(m,n));
          
          xhat = sparserandn(n,sp);  % true solution
          b_exact = A*xhat;          % exact data
          b = b_exact;
          
          tol_resAbz = 1e-7;
          tol_resATz = 1e-7; 

            
     
     case 'rank deficient, consistent, no noise'      

          rank = min(m,n)/2;
          sing_values_data.distribution = 'normal';
          sing_values_data.mu = 0; sing_values_data.stddev = 1;
          A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
          rel_cond_A = compute_rel_cond(A,min(m,n));
          
          xhat = sparserandn(n,sp);  % true solution
          b_exact = A*xhat;          % exact data
          b = b_exact;
          
          tol_resAbz = 1e-7;
          tol_resATz = 1e-7; 
          
          
          
        
     case 'rank-deficient, medium noise in R(A) complement'
       
          % m<=n and least-squares solution shall be not unique (no full rank)          
          
          noise_factor_rangeA_ortho = 0.5;
          rank = round(min(m,n)/2);   % must be smaller than m = min(m,n)
          real_setting = true;
          sing_values_data.distribution = 'normal';
          sing_values_data.mu = 0; sing_values_data.stddev = 1;
          
          A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
                 
          xhat = sparserandn(n,sp);  % true solution          
          b_exact = A*xhat;               % exact data
          noiselev = noise_factor_rangeA_ortho *norm(b_exact);
          
          % create offset with noise in the complement of R(A)
          N = null(A');             % null(A') = orthogonal complement of R(A)
          v = randn(size(N,2),1);
          noise_rangeA_ortho = N*v; 
          noise_rangeA_ortho = noiselev *noise_rangeA_ortho /norm(noise_rangeA_ortho);
          b = b_exact + noise_rangeA_ortho;    
          
          tol_resAbz = 1e-4;
          tol_resATz = 1e-6;

          
   
          
          
     case 'rank-deficient, large noise in R(A) complement'
       
          % m<=n and least-squares solution shall be not unique (no full rank)          
          
          noise_factor_rangeA_ortho = 5;
          rank = round(min(m,n)/2);   % must be smaller than m = min(m,n)
          real_setting = true;
          sing_values_data.distribution = 'normal';
          sing_values_data.mu = 0; sing_values_data.stddev = 1;
          
          A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
                 
          xhat = sparserandn(n,sp);  % true solution          
          b_exact = A*xhat;               % exact data
          noiselev = noise_factor_rangeA_ortho *norm(b_exact);
          
          % create offset with noise in the complement of R(A)
          N = null(A');             % null(A') = orthogonal complement of R(A)
          v = randn(size(N,2),1);
          noise_rangeA_ortho = N*v; 
          noise_rangeA_ortho = noiselev *noise_rangeA_ortho /norm(noise_rangeA_ortho);
          b = b_exact + noise_rangeA_ortho;    
          
          tol_resAbz = 1e-3;
          tol_resATz = 1e-5;
          
          
   
          


          
      case  'rank deficient, noise split into R(A) and R(A) complement'
  
          % m<=n and least-squares solution shall be not unique (no full rank), finally also add noise to b         
          
          noise_factor_rangeA_ortho = 1;
          noise_factor_rangeA = 0.1;
          rank = round(min(m,n)/2);   % must be smaller than m = min(m,n)
          
          real_setting = true;
          sing_values_data.distribution = 'normal';
          sing_values_data.mu = 0; sing_values_data.stddev = 1;
          
          A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
                  
          xhat = sparserandn(n,sp);  % true solution          
          b_exact = A*xhat;               % exact data
          noiselev_rangeA_ortho = noise_factor_rangeA_ortho *norm(b_exact);
          noiselev_rangeA = noise_factor_rangeA *norm(b_exact);
          
          % create noisy b 

          N = null(A');             % null(A') = orthogonal complement of R(A)
          v = randn(size(N,2),1);
          rangeA_ortho = N*v;
          rangeA_ortho = noiselev_rangeA_ortho* rangeA_ortho/norm(rangeA_ortho);
          
          noise_rangeA = A* randn(n,1);
          noise_rangeA = noiselev_rangeA/ norm(noise_rangeA) *noise_rangeA;

          b = b_exact + rangeA_ortho + noise_rangeA;      

          tol_resAbz = noiselev_rangeA;  % 1e-3
          tol_resATz = noiselev_rangeA;  % 1e-(6.5)
          
          
          
          
    case 'rank deficient, only noise in R(A) complement, well conditioned A'
      
      noiselev_rangeA_ortho = 0.5;
      rank = round(min(m,n)/2);
      
      sing_values_data.distribution = 'uniform';
      sing_values_data.sigma_min = 1; sing_values_data.sigma_max = 100;

      A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
      rel_cond_A = compute_rel_cond(A,rank);
          
      xhat = sparserandn(n,sp);  % true solution
      b_exact = A*xhat;               % exact data
      N = null(A');             % null(A') = orthogonal complement of R(A)
      v = randn(size(N,2),1);
      noise_rangeA_ortho = N*v;
      noise_rangeA_ortho = noiselev_rangeA_ortho* noise_rangeA_ortho/norm(noise_rangeA_ortho);
      b = b_exact + noise_rangeA_ortho;

      tol_resAbz = 1e-14;
      tol_resATz = 1e-14;
          
 
    
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

      tol_resAbz = 1e-12;
      tol_resATz = 1e-12;
          
          
      
      
      
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
      xhat = random_xhat_uniform_from_sp_min_max(sp, xhat_min, xhat_max, real_setting);

        
      % create noisy b 
      
      b_exact = A*xhat;               % exact data
      N = null(A');             % null(A') = orthogonal complement of R(A)
      v = randn(size(N,2),1);
      rangeA_ortho = N*v;
      rangeA_ortho = noiselev_rangeA_ortho* rangeA_ortho/norm(rangeA_ortho);
      noise_rangeA = A* randn(n,1);
      noise_rangeA = noiselev_rangeA/ norm(noise_rangeA) *noise_rangeA;
      
      b = b_exact + rangeA_ortho + noise_rangeA;      
          
      tol_resAbz = noiselev_rangeA;
      tol_resATz = 0.1;
      
      
      
      
     case 'dctmatrix, noise split into R(A) and R(A) complement'
         
       noiselev_rangeA_ortho = 5;
       noiselev_rangeA = 0.1;  % auch mal 0.1
       rank = min(m,n)/2;
         
       U=randn(m,rank); V=[eye(rank);dctmtx(rank)];
       A=U*V';
       rel_cond_A = compute_rel_cond(A,rank);
       
       xhat_min = 1;
       xhat_max = 6;
       xhat = random_xhat_uniform_from_sp_min_max(sp, xhat_min, xhat_max, real_setting);
       
       % create noisy b 
      
       b_exact = A*xhat;         % exact data
       N = null(A');             % null(A') = orthogonal complement of R(A)
       v = randn(size(N,2),1);
       rangeA_ortho = N*v;
       rangeA_ortho = noiselev_rangeA_ortho* rangeA_ortho/norm(rangeA_ortho);
       noise_rangeA = A* randn(n,1);
       noise_rangeA = noiselev_rangeA/ norm(noise_rangeA) *noise_rangeA;
      
       b = b_exact + rangeA_ortho + noise_rangeA;      
          
       tol_resAbz = noiselev_rangeA;
       tol_resATz = 0.1;
      
       
       
       
      % TODO
      case 'dctmatrix, impulsive noise'
       
       noiselev_rangeA_ortho = 5;
       noiselev_rangeA = 0.1;  % auch mal 0.1
       rank = min(m,n)/2;
         
       U=randn(m,rank); V=[eye(rank);dctmtx(rank)];
       A=U*V';
       rel_cond_A = compute_rel_cond(A,rank);
       
       xhat_min = 1;
       xhat_max = 6;
       xhat = random_xhat_uniform_from_sp_min_max(sp, xhat_min, xhat_max, real_setting);
       
       
       
       
      case 'full rank, small gaussian noise'
   
          noiselev = 0.1;
     
          A = randn(m,n);           % Operator
          rel_cond_A = compute_rel_cond(A,min(m,n));
          xhat = sparserandn(n,sp);  % true solution
          b_exact = A*xhat;               % exact data
          
          noise = randn(size(b_exact));
          b = b_exact + noiselev*noise/norm(noise);
          
          tol_resAbz = 1e-7;
          tol_resATz = 0.1;  
          
          
          
    case 'full rank, medium gaussian noise'
      
          noiselev = 0.5;
          
          A = randn(m,n);           % Operator
          rel_cond_A = compute_rel_cond(A,min(m,n));
          xhat = sparserandn(n,sp);  % true solution
          b_exact = A*xhat;               % exact data
          
          noise = randn(size(b_exact));
          b = b_exact + noiselev*noise/norm(noise);

          tol_resAbz = 1e-7;
          tol_resATz = 0.5;  
          
          
          
     case 'full rank, large gaussian noise'
      
          noiselev = 2;
          
          A = randn(m,n);           % Operator
          rel_cond_A = compute_rel_cond(A,min(m,n));
          xhat = sparserandn(n,sp);  % true solution
          b_exact = A*xhat;               % exact data
          
          noise = randn(size(b_exact));
          b = b_exact + noiselev*noise/norm(noise);
          
          tol_resAbz = 1e-7;
          tol_resATz = 0.5;            
          
          
          
    case 'large uniform noise, full rank'
      
          noiselev = 0.5;
          
          A = randn(m,n);           % Operator
          rel_cond_A = compute_rel_cond(A,min(m,n));
          xhat = sparserandn(n,sp);  % true solution
          b_exact = A*xhat;               % exact data
          
          noise = noiselev* (1 - 2*rand(size(b_exact)));
          b = b_exact + noise;

          tol_resAbz = 1e-7;
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

          tol_resAbz = 1e-7;
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
         
          tol_resAbz = 1e-7;
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

          tol_resAbz = 1e-7;
          tol_resATz = 0.5; 
          
          
          
  case 'rank-deficient, large noise in R(A) complement, complex'
          
          % m<=n and least-squares solution shall be not unique (no full rank)          
          
          noise_factor_rangeA_ortho = 0.5;
          rank = round(min(m,n)/2);   % must be smaller than m = min(m,n)
          real_setting = false;
          sing_values_data.distribution = 'normal';
          sing_values_data.mu = 0; sing_values_data.stddev = 1;
          
          A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
                 
          xhat = sparserandn(n,sp);  % true solution          
          b_exact = A*xhat;               % exact data
          noiselev = noise_factor_rangeA_ortho *norm(b_exact);
          
          % create offset with noise in the complement of R(A)
          N = null(A');             % null(A') = orthogonal complement of R(A)
          v = randn(size(N,2),1);
          noise_rangeA_ortho = N*v; 
          noise_rangeA_ortho = noiselev *noise_rangeA_ortho /norm(noise_rangeA_ortho);
          b = b_exact + noise_rangeA_ortho;    
          
          tol_resAbz = 1e-4;
          tol_resATz = 1e-6;
    
          
          
   case  'Franks example, complex_impulsive_noise_and_add_noise'
          
          num_comp_noise = ceil(m/40);
          noiselev_impulsive_noise_factor = 40;
          noiselev_add_noise_factor = 0.01;

          rank = 300;          
          real_setting = false;
          sing_values_data.distribution = 'uniform';
          sing_values_data.sigma_min = 1; sing_values_data.sigma_max = 2; 
          A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);

          xhat_min = 1;
          xhat_max = 6;
          xhat = random_xhat_uniform_from_sp_min_max(sp, xhat_min, xhat_max, real_setting);
          
          b_exact = A*xhat;
          b_max = noiselev_impulsive_noise_factor *max(abs(b_exact));
          supp_impulsive_noise = randi([1,m],1,num_comp_noise);
          
          b = b_exact;
          b(supp_impulsive_noise) = b(supp_impulsive_noise) + b_max*((rand(length(supp_impulsive_noise),1)-0.5)...
                                                            +1i*(rand(length(supp_impulsive_noise),1)-0.5));
          
          delta = noiselev_add_noise_factor *norm(b_exact);
          b_add = randn(size(b_exact));
          b = b + delta *b_add/norm(b_add);

          rel_cond_A = compute_rel_cond(A,rank);
          
          tol_resAbz = delta; 
          tol_resATz = 0.1;
          
          
          
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
 
 
          
  case 'impulsive noise, rank deficient, normalsqr, well-conditioned'
  
          num_comp_noise = ceil(min(m,n)/20);
          noiselev_impulsive = 10;
          
          rank = round(min(m,n)/2);
          sing_values_data.distribution = 'normalsqr';
          sing_values_data.mu = 0;
          sing_values_data.stddev = 0.5;
          

          A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
          xhat = sparserandn(n,sp);  % true solution          
          b_exact = A*xhat;               % exact data
          
          num_comp_noise = min(m,num_comp_noise);
          perm = randperm(m);
          b = b_exact;
          
          noise = randn(num_comp_noise,1);
          b(perm(1:num_comp_noise)) = b(perm(1:num_comp_noise)) + noiselev_impulsive* noise / norm(noise); 
          
          tol_resAbz = noiselev_impulsive*1e-6;
          tol_resATz = noiselev_impulsive*1e-3; 

          
          
          
  case 'impulsive noise, rank deficient, uniform, well-conditioned'
  
          num_comp_noise = ceil(min(m,n)/20);
          noiselev_impulsive = 10;
          
          rank = round(min(m,n)/2);
          sing_values_data.distribution = 'uniform';
          sing_values_data.sigma_min = 1e-3;
          sing_values_data.sigma_max = 3;
          

          A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
          xhat = sparserandn(n,sp);  % true solution          
          b_exact = A*xhat;               % exact data
          
          num_comp_noise = min(m,num_comp_noise);
          perm = randperm(m);
          b = b_exact;
          
          noise = randn(num_comp_noise,1);
          b(perm(1:num_comp_noise)) = b(perm(1:num_comp_noise)) + noiselev_impulsive* noise / norm(noise); 
          
          tol_resAbz = noiselev_impulsive*1e-6;
          tol_resATz = noiselev_impulsive*1e-3; 
          
          
    
          
  case 'impulsive noise, rank-deficient, well conditioned A'
  
          num_comp_noise = 10;
          noiselev_impulsive_noise_factor = 0.1;
          
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
          
          tol_resAbz = noiselev_impulsive_noise_factor*noiselev_impulsive*1e-5;
          tol_resATz = noiselev_impulsive_noise_factor*noiselev_impulsive*1e-3; 
          
          
          
          
    case 'impulsive noise, rank-deficient, A with unif distr sv, well-conditioned A and xhat'
  
          num_comp_noise = ceil(min(m,n)/20);
          noiselev_impulsive_noise_factor = 10;

          rank = round(min(m,n)/2);
          sing_values_data.distribution = 'uniform';
          sing_values_data.sigma_min = 1;
          sing_values_data.sigma_max = 2;
          
          xhat_min = 0.1;
          xhat_max = 1;

          A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
                    
          xhat = random_xhat_uniform_from_sp_min_max(sp, xhat_min, xhat_max, real_setting);        
          b_exact = A*xhat;               % exact data
          
          num_comp_noise = min(m,num_comp_noise);
          perm = randperm(m);
          b = b_exact;
          
          noiselev_impulsive = noiselev_impulsive_noise_factor *norm(b);
          noise = randn(num_comp_noise,1);
          b(perm(1:num_comp_noise)) = b(perm(1:num_comp_noise)) + noiselev_impulsive* noise / norm(noise); 
          
          tol_resAbz = 1e-4 *noiselev_impulsive;
          tol_resATz = 1e-5;         

          
          
          
     case 'impulsive noise, rank-deficient, A with unif distr sv, medium-conditioned A'
  
          num_comp_noise = min(m,n)/20;
          noiselev_impulsive_noise_factor = 10;

          rank = round(min(m,n)/2);
          sing_values_data.distribution = 'uniform';
          sing_values_data.sigma_min = 0.01;
          sing_values_data.sigma_max = 1;
          
          xhat_min = 0.1;
          xhat_max = 1;

          A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
                    
          xhat = random_xhat_uniform_from_sp_min_max(sp, xhat_min, xhat_max, real_setting);        
          b_exact = A*xhat;               % exact data
          
          num_comp_noise = min(m,num_comp_noise);
          perm = randperm(m);
          b = b_exact;
          
          noiselev_impulsive = noiselev_impulsive_noise_factor *norm(b);
          noise = randn(num_comp_noise,1);
          b(perm(1:num_comp_noise)) = b(perm(1:num_comp_noise)) + noiselev_impulsive* noise / norm(noise); 
          
          tol_resAbz = 1e-3*noiselev_impulsive;
          tol_resATz = 0.5*1e-4*noiselev_impulsive;  
          
          

     case 'impulsive noise, rank-deficient, A with unif distr sv, bad-conditioned A'
  
          num_comp_noise = min(m,n)/20;
          noiselev_impulsive_noise_factor = 10;

          rank = round(min(m,n)/2);
          sing_values_data.distribution = 'uniform';
          sing_values_data.sigma_min = 0.001;
          sing_values_data.sigma_max = 1;
          
          xhat_min = 0.1;
          xhat_max = 1;

          A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
                    
          xhat = random_xhat_uniform_from_sp_min_max(sp, xhat_min, xhat_max, real_setting);        
          b_exact = A*xhat;               % exact data
          
          num_comp_noise = min(m,num_comp_noise);
          perm = randperm(m);
          b = b_exact;
          
          noiselev_impulsive = noiselev_impulsive_noise_factor *norm(b);
          noise = randn(num_comp_noise,1);
          b(perm(1:num_comp_noise)) = b(perm(1:num_comp_noise)) + noiselev_impulsive* noise / norm(noise); 
          
          tol_resAbz = 1e-3*noiselev_impulsive;
          tol_resATz = 0.5*1e-4*noiselev_impulsive;  

          
          
     case 'less impulsive noise, rank-deficient, A with unif distr sv, bad-conditioned A'
  
          num_comp_noise = min(m,n)/20;
          noiselev_impulsive_noise_factor = 2;

          rank = round(min(m,n)/2);
          sing_values_data.distribution = 'uniform';
          sing_values_data.sigma_min = 0.001;
          sing_values_data.sigma_max = 1;
          
          xhat_min = 0.1;
          xhat_max = 1;

          A = random_rank_deficient_matrix_with_condition(m,n,rank,real_setting,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
                    
          xhat = random_xhat_uniform_from_sp_min_max(sp, xhat_min, xhat_max, real_setting);        
          b_exact = A*xhat;               % exact data
          
          num_comp_noise = min(m,num_comp_noise);
          perm = randperm(m);
          b = b_exact;
          
          noiselev_impulsive = noiselev_impulsive_noise_factor *norm(b);
          noise = randn(num_comp_noise,1);
          b(perm(1:num_comp_noise)) = b(perm(1:num_comp_noise)) + noiselev_impulsive* noise / norm(noise); 
          
          tol_resAbz = 1e-3*noiselev_impulsive;
          tol_resATz = 0.5*1e-4*noiselev_impulsive; 
          

   case 'impulsive noise and 10% additional noise, rank-deficient, bad-conditioned A and xhat'
  
          num_comp_noise = 10;
          noiselev_impulsive_noise_factor = 10;
          noiselev_add_noise_factor = 0.1;

          rank = round(min(m,n)/2);
          sing_values_data.distribution = 'uniform';
          sing_values_data.sigma_min = 0.001;
          sing_values_data.sigma_max = 1;

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
          
          add_noise = randn(size(b));
          b = b + noiselev_add_noise_factor * add_noise/norm(add_noise) * norm(b);
          
          tol_resAbz = 1e-5*noiselev_impulsive_noise_factor *noiselev_impulsive;
          tol_resATz = 1e-4*noiselev_impulsive_noise_factor *noiselev_impulsive;              
 
          
          
     case 'impulsive noise and 1% additional noise, rank-deficient, bad-conditioned A and xhat'
  
          num_comp_noise = 10;
          noiselev_impulsive_noise_factor = 10;
          noiselev_add_noise_factor = 0.01;

          rank = round(min(m,n)/2);
          sing_values_data.distribution = 'uniform';
          sing_values_data.sigma_min = 0.001;
          sing_values_data.sigma_max = 1;

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
          
          add_noise = randn(size(b));
          b = b + noiselev_add_noise_factor * add_noise/norm(add_noise) * norm(b);
          
          tol_resAbz = 1e-5*noiselev_impulsive_noise_factor *noiselev_impulsive;
          tol_resATz = 1e-7*noiselev_impulsive_noise_factor *noiselev_impulsive;              
 
          

    
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
         
          tol_resAbz = noiselev_gaussian_noise_factor;
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
      elseif strcmp(sing_values_data.distribution,'normalsqr')
          sing_values = sing_values_data.mu + sing_values_data.stddev *randn(1,rank);
          sing_values = sing_values.^2;
      else
          error('Failed to create system matrix: No matching string!')
      end
      
      D = diag(sing_values);  
      A=U*D*V';
      
  end
      
      
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

  function xhat = random_xhat_uniform_from_sp_min_max(sp, xhat_min, xhat_max, real_setting)
      xhat=zeros(n,1);
      rpn = randperm(n); 
      supp = rpn(1:sp);
      if real_setting
          xhat(supp)= xhat_min+ (xhat_max-xhat_min)*rand(sp,1);
      else
          xhat(supp)= xhat_min+ (xhat_max-xhat_min)*(rand(sp,1) +1i*rand(sp,1));
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
