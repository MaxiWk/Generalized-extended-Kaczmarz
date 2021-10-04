function problem_data = set_up_instance(m,n,sp,real_setting,experiment_description)
  
  b_exact = 0;  
  tol_resA = 0;
  tol_resAT = 0;
  
  switch experiment_description
    
    
    
  case 'small_example'
    
          A = [1,2; 1,0; 0,1];
          b = [0;1;1];
          xhat = [0.5;0];
    
     case 'solvable, no noise'
       
          A = randn(m,n);           % Operator
          xhat = sparserandn(n,sp);  % true solution
          b = A*xhat;               % exact data


          
     case 'not solvable, noise only in the complement of range'
       
          noiselev_rangeA_ortho = 0.5;
          
          A = randn(m,n);           % Operator
          xhat = sparserandn(n,sp);  % true solution
          b_exact = A*xhat;               % exact data
          N = null(A');             % null(A') = orthogonal complement of R(A)
          v = randn(size(N,2),1);
          noise_rangeA_ortho = N*v;
          noise_rangeA_ortho = noiselev_rangeA_ortho* noise_rangeA_ortho/norm(noise_rangeA_ortho);
          b = b_exact + noise_rangeA_ortho;
          
          
          
    case 'not solvable, noise only in the complement of range, well conditioned A'
       
      noiselev_rangeA_ortho = 0.5;
      tol_resA = 1e-3;
      tol_resAT = 1e-3;
      sigma_min = 1;
      sigma_max = 2;
      rank = round(2*min(m,n)/3);

      if real_setting
          U=orth(randn(m,rank)); V=orth(randn(n,rank));
      else
          U=orth(randn(m,rank)+1i*randn(m,r)); V=orth(randn(n,rank)+1i*randn(n,rank));
      end
      D=diag( sigma_min + (sigma_max-sigma_min)*rand(1,rank) );
      A=U*D*V';
          
      xhat = sparserandn(n,sp);  % true solution
      b_exact = A*xhat;               % exact data
      N = null(A');             % null(A') = orthogonal complement of R(A)
      v = randn(size(N,2),1);
      noise_rangeA_ortho = N*v;
      noise_rangeA_ortho = noiselev_rangeA_ortho* noise_rangeA_ortho/norm(noise_rangeA_ortho);
      b = b_exact + noise_rangeA_ortho;

          
 
    
    case 'not solvable, noise in the complement of range, well conditioned A and xhat'  
        
      noiselev_rangeA_ortho = 5;
      tol_resA = 1e-5;
      tol_resAT = 1e-5;
      sigma_min = 1;
      sigma_max = 2;
      xhat_min = 1;
      xhat_max = 6;
      rank = round(2*min(m,n)/3);

      if real_setting
          U=orth(randn(m,rank)); V=orth(randn(n,rank));
      else
          U=orth(randn(m,rank)+1i*randn(m,r)); V=orth(randn(n,rank)+1i*randn(n,rank));
      end
      D=diag( sigma_min + (sigma_max-sigma_min)*rand(1,rank) );
      A=U*D*V';
                          
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

      
      
      
    case 'not solvable, noise in the complement of range and noise in range, well conditioned A and xhat'  
        
      noiselev_rangeA_ortho = 5;
      noiselev_rangeA = 0.01;
      sigma_min = 1;
      sigma_max = 2;
      xhat_min = 1;
      xhat_max = 6;
      rank = round(2*min(m,n)/3);

      if real_setting
          U=orth(randn(m,rank)); V=orth(randn(n,rank));
      else
          U=orth(randn(m,rank)+1i*randn(m,r)); V=orth(randn(n,rank)+1i*randn(n,rank));
      end
      D=diag( sigma_min + (sigma_max-sigma_min)*rand(1,rank) );
      A=U*D*V';
                          
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
      noise_rangeA = A* randn(n,1);
      noise_rangeA = noiselev_rangeA/ norm(noise_rangeA) *noise_rangeA;
      
      b = b_exact + rangeA_ortho + noise_rangeA;      
          
      tol_resA = noiselev_rangeA *norm(b_exact);
      tol_resAT = 0.1;
      
          
      case 'not solvable, little gaussian noise'
   
          noiselev = 0.1;
     
          A = randn(m,n);           % Operator
          xhat = sparserandn(n,sp);  % true solution
          b_exact = A*xhat;               % exact data
          
          noise = randn(size(b));
          b = b_exact + noiselev*noise/norm(noise);
          
    case 'not solvable, medium gaussian noise'
      
          noiselev = 0.5;
          
          A = randn(m,n);           % Operator
          xhat = sparserandn(n,sp);  % true solution
          b_exact = A*xhat;               % exact data
          
          noise = randn(size(b_exact));
          b = b_exact + noiselev*noise/norm(noise);
          
     case 'not solvable, gaussian noise'
      
          noiselev = 2;
          
          A = randn(m,n);           % Operator
          xhat = sparserandn(n,sp);  % true solution
          b_exact = A*xhat;               % exact data
          
          noise = randn(size(b_exact));
          b = b_exact + noiselev*noise/norm(noise);
          
    case 'large uniform noise'
      
          noiselev = 0.5;
          
          A = randn(m,n);           % Operator
          xhat = sparserandn(n,sp);  % true solution
          b_exact = A*xhat;               % exact data
          
          noise = noiselev* (1 - 2*rand(size(b_exact)));
          b = b_exact + noise;




    case 'rank deficient, small uniform noise'
      
          noiselev = 0.01;
          rank = round(2*min(m,n)/3); 
          A = random_rank_deficient_matrix(m,n,rank);
          xhat = sparserandn(n,sp);
          b_exact = A*xhat;
          
          noise = noiselev* (1 - 2*rand(size(b_exact)));
          b = b_exact + noise;         

 




    case 'rank deficient, small gaussian noise'
      
          noiselev = 0.01;
          rank = round(2*min(m,n)/3); 
          A = random_rank_deficient_matrix(m,n,rank);
          xhat = sparserandn(n,sp);
          b_exact = A*xhat;
          
          noise = noiselev* randn(size(b_exact));
          b = b_exact + noise;   
         
         
          
          
        
    case 'rank deficient, medium uniform noise'
    
          noiselev = 0.1;
          rank = round(2*min(m,n)/3);
          A = random_rank_deficient_matrix(m,n,rank);
          xhat = sparserandn(n,sp);
          b_exact = A*xhat;
          
          noise = noiselev* (1 - 2*rand(size(b_exact)));
          b = b_exact + noise;





    case 'rank deficient, large b plus medium uniform noise'
    
          noiselev = 0.1;
          rank = round(2*min(m,n)/3);
          A = random_rank_deficient_matrix(m,n,rank);
          xhat = sparserandn(n,sp);  
          xhat = xhat / norm(xhat) * 10;
          b_exact = A*xhat;
          
          noise = noiselev* (1 - 2*rand(size(b_exact)));
          b = b_exact + noise;          
          
          
          
          
    case 'not solvable, noise in complement of range and gaussian noise added'
       
          noiselev_rangeA_ortho = 0.5;
          noiselev = 0.5;
          
          A = randn(m,n);           % Operator
          xhat = sparserandn(n,sp);  % true solution
          b = A*xhat;               % exact data
          N = null(A');             % null(A') = orthogonal complement of R(A)
          v = randn(size(N,2),1);
          noise_rangeA_ortho = N*v;
          noise_rangeA_ortho = noiselev_rangeA_ortho* noise_rangeA_ortho/norm(noise_rangeA_ortho);
          b = b + noise_rangeA_ortho;
          b = b + noiselev*randn(m,1);
          
          
    case 'not solvable, noise in complement of range and uniform noise added'
       
          noiselev_rangeA_ortho = 0.5;
          noiselev = 0.5;
          
          A = randn(m,n);           % Operator
          xhat = sparserandn(n,sp);  % true solution
          b = A*xhat;               % exact data
          N = null(A');             % null(A') = orthogonal complement of R(A)
          v = randn(size(N,2),1);
          noise_rangeA_ortho = N*v;
          noise_rangeA_ortho = noiselev_rangeA_ortho* noise_rangeA_ortho/norm(noise_rangeA_ortho);
          b = b + noise_rangeA_ortho;
          b = b + noiselev*rand(m,1);
          
     
    
 case 'rank-deficient, noise in complement of range'
       
          % m<=n and least-squares solution shall be not unique (no full rank)          
          
          noiselev_rangeA_ortho = 0.5;
          rank = round(min(m,n)/2);   % must be smaller than m = min(m,n)
          
          A = random_rank_deficient_matrix(m,n,rank);
                  
          xhat = sparserandn(n,sp);  % true solution          
          b_exact = A*xhat;               % exact data
          
          % create offset with noise in the complement of R(A)
          N = null(A');             % null(A') = orthogonal complement of R(A)
          v = randn(size(N,2),1);
          noise_rangeA_ortho = N*v; 
          noise_rangeA_ortho = noiselev_rangeA_ortho* noise_rangeA_ortho/norm(noise_rangeA_ortho);
          b = b_exact + noise_rangeA_ortho;    
          
          % then: pinv(A)*b = pinv(A)*b_exact 
          
          
 

  case  'rank-deficient, gaussian noise in complement of range and then uniform noise added' 
  
          % m<=n and least-squares solution shall be not unique (no full rank), finally also add noise to b         
          
          noiselev_rangeA_ortho = 0.5;
          add_noise_norm = 0.1;
          rank = round(min(m,n)/2);   % must be smaller than m = min(m,n)
          
          A = random_rank_deficient_matrix(m,n,rank);
                  
          xhat = sparserandn(n,sp);  % true solution          
          b_exact = A*xhat;               % exact data
          
          % create offset with noise in the complement of R(A)
          N = null(A');             % null(A') = orthogonal complement of R(A)
          v = randn(size(N,2),1);
          noise_rangeA_ortho = N*v; 
          noise_rangeA_ortho = noiselev_rangeA_ortho* noise_rangeA_ortho/norm(noise_rangeA_ortho);
          b = b_exact + noise_rangeA_ortho; 
          
          b = b + add_noise_norm*rand(m,1);
          
          
         
         
          
          
  case 'large impulsive noise'
  
          num_comp_noise = 10;

          %A = randn(m,n);
          rank = floor(min(m,n)/2); 
          A = randn(m,n);                  
          xhat = sparserandn(n,sp);  % true solution          
          b_exact = A*xhat;               % exact data
          
          num_comp_noise = min(m,num_comp_noise);
          perm = randperm(m);
          b = b_exact;
          noise = randn(num_comp_noise,1);
          noise_corr_factor = 0.1* norm(b)/norm(noise);
          b(perm(1:num_comp_noise)) = b(perm(1:num_comp_noise)) + noise_corr_factor* noise; 
          

          
          
          
          
          
     otherwise 
          disp('No experiment chosen!')
 
          
         
  end
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function A = random_rank_deficient_matrix(m,n,rank)
    
          swap_mn = false;      % wlog m<n (swap m and n twice)
          if n < m
            [m,n] = deal(n,m);   
            swap_mn = true;
          end
          
          % generate random orthogonal matrices U and V, A will be given by svd
          tmp = randn(m,m);
          [U,~] = qr(tmp);
          tmp = randn(n,n);
          [V,~] = qr(tmp);
          
          % build up A with a rank deficient singular value decomposition
          sing_values = sparserandn(m, rank); 
          Sigma = [diag(sing_values), zeros(m,n-m)];        
          A = U * Sigma * V';
          
          if swap_mn
            [m,n] = deal(n,m);   % swap m and n back
            A = A';
          end
          
  end
  



 problem_data.A = A; 
 problem_data.xhat = xhat;
 problem_data.b_exact = b_exact;
 problem_data.b = b;
 problem_data.tol_resA = tol_resA;
 problem_data.tol_resAT = tol_resAT;

end
