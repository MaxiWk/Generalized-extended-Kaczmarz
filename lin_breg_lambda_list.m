lambda_list = 2.567914930816152e-02 + [0, 0.01, 0.1,0.2,0.3,0.4,0.5,...
                                       2, 3, 4, 100] ;

lambda_list = 0.03 + [0, 0.01, 0.1,0.2,0.3,0.4,0.5,...
                                       2, 3, 4, 100] ;                                      
rand('state',0);
randn('state',0);


m = 2;
n = 5;


maxiter = 1e5; % Number of iterations 
time_per_repeat = 5; % time to run if 'experiment_time'
iter_save = ceil(maxiter/500);  % each such number of iterations, a data point is added in the error plot

rowsamp = 'rownorms squared';
colsamp = 'colnorms squared'; 



A = randn(m,n);
b = rand(m,1);

xplus = pinv(A)*b;
disp('Moore penrose inverse:')
disp(xplus)



x_list_lin_bregman = zeros(length(xplus),length(lambda_list));

ATA = A'*A;
ATb = A'*b;
lin_breg_stepsize = 1/norm(A)^2;




i = 1;

for lambda = lambda_list
  
  disp(['Linearized Bregman method, ' num2str(i) '/' num2str(length(lambda_list))])

  S = @(z) max(abs(z)-lambda,0).*sign(z); % Soft thresholding 

  xdual = zeros(n,1);
  x = zeros(n,1);

  for iter = 1:maxiter
    xdual = xdual - lin_breg_stepsize*(ATA*x-ATb);
    x = S(xdual);
  end

  x_list_lin_bregman(:,i) = x;
  
  i = i+1;
  
end