%{
A = [1,2,3; 1,0,2; -1,-2,-3];
xhat = [0;1;2];
b = A*xhat;

n = size(A,2);
%}

m = 200;
n = 500;
sp = 5;
real_setting = true;
problem_data = set_up_instance(m,n,sp,real_setting,'full rank, consistent, no noise');
A = problem_data.A;
b = problem_data.A;


% lin bregman with large lambda for l1 norm solution

lambda = 10;
maxiter = 1000;
S = @(z) max(abs(z)-lambda,0).*sign(z); % Soft thresholding

lin_breg_stepsize = 1/norm(A)^2;
x = zeros(n,1);
xdual = zeros(n,1);

for i = 1:maxiter
    xdual = xdual - lin_breg_stepsize*A'*(A*x-b);
    x = S(xdual);
    %disp(['Residual: ' num2str(norm(A*x-b))])
end

disp('Solution found with lin bregman method: ')
%disp(x)

disp(['Residual: ' num2str(norm(A*x-b))])
disp('')



% try out acceleration of lin bregman

lin_breg_stepsize = 1/norm(A)^2;
x = zeros(n,1);
xstar = zeros(n,1);
xtilde = zeros(n,1);

alpha = 0.5;
omega = 1;

for i = 1:maxiter
    xstar_old = xstar;
    xstar = xtilde - lin_breg_stepsize*A'*(A*x-b);
    xtilde = alpha*xstar + (1-alpha)*xstar_old;
    x = S(xtilde);
    %disp(['Residual: ' num2str(norm(A*x-b))])
end

disp('Solution found by accelerated lin bregman method: ')
%disp(x)

disp(['Residual: ' num2str(norm(A*x-b))])




%% compare sparse Kaczmarz and 'acceleration' idea

%%%%%%%% paramters and helper functions %%%%%%%%%%%%%

lambda = 1;
maxiter = 100;

S = @(z) max(abs(z)-lambda,0).*sign(z); % Soft thresholding

norma = sum(A.^2,2);
p = norma./sum(norma); 
P = cumsum(p);
samp_row = @(k) nnz(rand>P)+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% sparse Kaczmarz

x = zeros(n,1);
xstar = zeros(n,1);

rand('state',0);

for iter = 1:maxiter
    r = samp_row(iter);
    a = A(r,:);
    xstar = xstar - ( a*x - b(r) )/norm(a)^2 *a';
    x = S(xstar);
end

disp(['Sparse Kaczmarz, Residual: ' num2str(norm(A*x-b))])



% try new acceleration of sparse Kaczmarz

x = zeros(n,1);
xstar = zeros(n,1);
xtilde = zeros(n,1);

alpha = 0.1;
omega = 1;

rand('state',0);

for iter = 1:maxiter
    
    r = samp_row(iter);
    a = A(r,:);
    xstar_old = xstar;
    xstar = xtilde - omega *( a*x - b(r) )/norm(a)^2 *a';
    xtilde = alpha*xstar + (1-alpha)*xstar_old;
    x = S(xtilde);
    
    if mod(iter,100)==0
        disp(['Residual: ' num2str(norm(A*x-b))])
    end
    
end

