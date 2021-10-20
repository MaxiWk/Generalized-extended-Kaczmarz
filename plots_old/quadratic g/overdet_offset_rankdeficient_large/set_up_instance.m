function [A,b,b_exact,xhat] = set_up_instance(m,n,sp,experiment_description)
  
  b_exact = 0;  % just to not have it undefined
  
  switch experiment_description
    
     case 'solvable, no noise'
       
          A = randn(m,n);           % Operator
          xhat = sparserandn(n,sp);  % true solution
          b = A*xhat;               % exact data


          
     case 'not solvable, noise only in the complement of range'
       
          offset_norm = 0.5;
          
          A = randn(m,n);           % Operator
          xhat = sparserandn(n,sp);  % true solution
          b = A*xhat;               % exact data
          N = null(A');             % null(A') = orthogonal complement of R(A)
          v = randn(size(N,2),1);
          offset = N*v;
          offset = offset_norm* offset/norm(offset);
          b = b + offset;
 
          
          
      case 'not solvable, little gaussian noise'
   
          noiselev = 0.1;
     
          A = randn(m,n);           % Operator
          xhat = sparserandn(n,sp);  % true solution
          b = A*xhat;               % exact data
          
          noise = randn(size(b));
          b = b + noiselev*noise/norm(noise);
          
    case 'not solvable, medium gaussian noise'
      
          noiselev = 0.5;
          
          A = randn(m,n);           % Operator
          xhat = sparserandn(n,sp);  % true solution
          b = A*xhat;               % exact data
          
          noise = randn(size(b));
          b = b + noiselev*noise/norm(noise);
          
     case 'not solvable, gaussian noise'
      
          noiselev = 2;
          
          A = randn(m,n);           % Operator
          xhat = sparserandn(n,sp);  % true solution
          b = A*xhat;               % exact data
          
          noise = randn(size(b));
          b = b + noiselev*noise/norm(noise);
          
    case 'medium uniform noise'
      
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
      
          noiselev = 0.1;
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
          
          
         
          
          
    case 'not solvable, noise in complement of range and gaussian noise added'
       
          offset_norm = 0.5;
          noiselev = 0.5;
          
          A = randn(m,n);           % Operator
          xhat = sparserandn(n,sp);  % true solution
          b = A*xhat;               % exact data
          N = null(A');             % null(A') = orthogonal complement of R(A)
          v = randn(size(N,2),1);
          offset = N*v;
          offset = offset_norm* offset/norm(offset);
          b = b + offset;
          b = b + noiselev*randn(m,1);
          
          
    case 'not solvable, noise in complement of range and uniform noise added'
       
          offset_norm = 0.5;
          noiselev = 0.5;
          
          A = randn(m,n);           % Operator
          xhat = sparserandn(n,sp);  % true solution
          b = A*xhat;               % exact data
          N = null(A');             % null(A') = orthogonal complement of R(A)
          v = randn(size(N,2),1);
          offset = N*v;
          offset = offset_norm* offset/norm(offset);
          b = b + offset;
          b = b + noiselev*rand(m,1);
          
     
    
 case 'rank-deficient, noise in complement of range'
       
          % m<=n and least-squares solution shall be not unique (no full rank)          
          
          offset_norm = 0.5;
          rank = round(min(m,n)/2);   % must be smaller than m = min(m,n)
          
          A = random_rank_deficient_matrix(m,n,rank);
                  
          xhat = sparserandn(n,sp);  % true solution          
          b_exact = A*xhat;               % exact data
          
          % create offset with noise in the complement of R(A)
          N = null(A');             % null(A') = orthogonal complement of R(A)
          v = randn(size(N,2),1);
          offset = N*v; 
          offset = offset_norm* offset/norm(offset);
          b = b_exact + offset;    
          
          % then: pinv(A)*b = pinv(A)*b_exact 
          
          
 

  case  'rank-deficient, gaussian noise in complement of range and then uniform noise added' 
  
          % m<=n and least-squares solution shall be not unique (no full rank), finally also add noise to b         
          
          offset_norm = 0.5;
          add_noise_norm = 0.1;
          rank = round(min(m,n)/2);   % must be smaller than m = min(m,n)
          
          A = random_rank_deficient_matrix(m,n,rank);
                  
          xhat = sparserandn(n,sp);  % true solution          
          b_exact = A*xhat;               % exact data
          
          % create offset with noise in the complement of R(A)
          N = null(A');             % null(A') = orthogonal complement of R(A)
          v = randn(size(N,2),1);
          offset = N*v; 
          offset = offset_norm* offset/norm(offset);
          b = b_exact + offset; 
          
          b = b + add_noise_norm*rand(m,1);
          
          
         
          
          
          
   case 'small example, solvable, nonunique'
     
          randn('state',0);
          A = randn(2,5);
          xhat = sparserandn(5,2);
          b = A*xhat;
          xplus = pinv(A)*b;  % least squares sol with min 2 norm

          
   case 'small example, offset'
     
          randn('state',0);
          A = randn(2,5);
          xhat = sparserandn(5,2);          
          b = A*xhat;
          N = null(A');             % null(A') = orthogonal complement of R(A)
          v = randn(size(N,2),1);
          offset = N*v;
          offset = offset_norm* offset/norm(offset);
          b = b + offset;
          xplus = pinv(A)*b;  % least squares sol with min 2 norm
          

          
     otherwise 
          disp('No experiment chosen!')
 
         
  end
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  % helper function
  
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
  
  
  
end
