% Here, random problem instances (A, b, xhat) are generated.

function problem_data = set_up_instance(m, n, sp, real, experiment_description)
  
  bhat = 0;  
  tol_resAbz = 0;
  tol_resATz = 0;
  
  switch experiment_description
    
    

      
      
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % experiment (i) (sparse least squares examples)
          
   
     case 'rank-deficient, large noise in R(A) complement'
       
          noise_factor_rangeA_ortho = 5;
          rank = round(min(m,n)/2);   % must be smaller than m = min(m,n)
          sing_values_data.distribution = 'uniform';
          sing_values_data.sigma_min = 1e-3; 
          sing_values_data.sigma_max = 100;
          
          A = random_rank_deficient_matrix_with_condition(m,n,rank,real,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
                 
          xhat = sparserandn(n,sp,real);  % true solution          
          bhat = A*xhat;               % exact data
          
          noiselev = noise_factor_rangeA_ortho *norm(bhat);
          b = add_noise_in_RAc(bhat, A, noiselev);
          
          tol_resAbz = 1e-4;
          tol_resATz = 1e-6;
          
          
    

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % experiment (ii) (examples with impulsive noise)

   
   case 'fixed impulsive noise, rank deficient, uniform, medium conditioned'
  
          num_comp_noise = ceil(min(m,n)/20);
          noiselev_impulsive_factor = 5;
          
          rank = round(min(m,n)/2);
          sing_values_data.distribution = 'uniform';
          sing_values_data.sigma_min = 0.001; 
          sing_values_data.sigma_max = 10;
          
          A = random_rank_deficient_matrix_with_condition(m,n,rank,real,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
          xhat = sparserandn(n,sp,real);  % true solution          
          bhat = A*xhat;               % exact data          
          
          noiselev_impulsive = noiselev_impulsive_factor * max(abs(bhat));
          b = add_fixed_impulsive_noise(bhat, num_comp_noise, noiselev_impulsive, real);
                  
          tol_resAbz = noiselev_impulsive * 1e-6;
          tol_resATz = noiselev_impulsive * 1e-3; 
 
 
          
          
   case 'fixed impulsive noise, rank deficient, uniform, bad conditioned'
  
          num_comp_noise = ceil(min(m,n)/20);
          noiselev_impulsive_factor = 5;
          
          rank = round(min(m,n)/2);
          sing_values_data.distribution = 'uniform';
          sing_values_data.sigma_min = 1e-3;
          sing_values_data.sigma_max = 100;
          

          A = random_rank_deficient_matrix_with_condition(m,n,rank,real,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
          xhat = sparserandn(n,sp,real);  % true solution          
          bhat = A*xhat;               % exact data          
          
          noiselev_impulsive = noiselev_impulsive_factor *max(abs(bhat));
          b = add_fixed_impulsive_noise(bhat, num_comp_noise, noiselev_impulsive, real);
                  
          tol_resAbz = noiselev_impulsive*1e-6;
          tol_resATz = noiselev_impulsive*1e-3; 
           
          
          
          
     case 'fixed small impulsive noise, rank deficient, uniform, medium conditioned'      
  
          num_comp_noise = ceil(min(m,n)/20);
          noiselev_impulsive_factor = 0.1;
          
          rank = round(min(m,n)/2);
          sing_values_data.distribution = 'uniform';
          sing_values_data.sigma_min = 0.001; 
          sing_values_data.sigma_max = 10;
          

          A = random_rank_deficient_matrix_with_condition(m,n,rank,real,sing_values_data);
          rel_cond_A = compute_rel_cond(A,rank);
          xhat = sparserandn(n,sp,real);  % true solution          
          bhat = A*xhat;               % exact data          
          
          noiselev_impulsive = noiselev_impulsive_factor *max(abs(bhat));
          b = add_fixed_impulsive_noise(bhat, num_comp_noise, noiselev_impulsive, real);
                  
          tol_resAbz = noiselev_impulsive*1e-6;
          tol_resATz = noiselev_impulsive*1e-3; 
          
          
          
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


  function b = add_noise_in_RAc(b,A,noiselev)
      ONB_A_adj = null(A');            
      noise_coeff = randn(size(ONB_A_adj,2),1);
      noise_in_RAc = ONB_A_adj*noise_coeff; 
      noise_in_RAc = noiselev *noise_in_RAc /norm(noise_in_RAc);
      b = b + noise_in_RAc;  
  end


  function b = add_fixed_impulsive_noise(b, num_comp_noise, noiselev, real_setting)
          
      supp_impulsive_noise = rand_sample(m, num_comp_noise);
      rand_signs = (-1).^ randi(2, 1, num_comp_noise)';
      sparse_random_vector = rand_signs * noiselev;

      if ~real_setting
          rand_signs = (-1).^ randi(2, 1, num_comp_noise)';
          sparse_random_vector = sparse_random_vector + 1i * rand_signs * noiselev;
          sparse_random_vector = sparse_random_vector / sqrt(2);
      end
      
      b(supp_impulsive_noise) = b(supp_impulsive_noise) + sparse_random_vector;

  end


  function b = add_randn_impulsive_noise(b, num_comp_noise, noiselev, real_setting)
          
      supp_impulsive_noise = rand_sample(m,num_comp_noise);
      sparse_random_vector = noiselev *(2*rand(length(supp_impulsive_noise),1)-1);

      if ~real_setting
          sparse_random_vector = sparse_random_vector + 1i *noiselev *(2*rand(length(supp_impulsive_noise),1)-1);
          sparse_random_vector = sparse_random_vector / sqrt(2);
      end

      b(supp_impulsive_noise) = b(supp_impulsive_noise) + sparse_random_vector;

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
 problem_data.bhat = bhat;
 problem_data.b = b;
 problem_data.tol_resAbz = tol_resAbz;
 problem_data.tol_resATz = tol_resATz;

end




