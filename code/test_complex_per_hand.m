%% build up 2x2 example 

%A = [1+2i,2-i;3i,4];
A = [1+2i, 2-i; 2+4i, 4-2i];
xhat = [5,1]';
b = A*xhat;


%% build up larger random example

m = 50;
n = 20;
A = randn(m,n) + 1i*randn(m,n);
xhat = randn(size(A,2),1);
b = A*xhat + randn(size(A,1),1);


%% define sampling functions

% define sampling functions

norma = sum(abs(A).^2,2);      % squared norms of rows of A

norma_col = sum(abs(A).^2,1);  % squared norms of cols of A

% sampling function for rows
p = norma./sum(norma); 
P = cumsum(p);
samp_row = @(k) nnz(rand>P)+1;

% sampling function for columns
pcol = norma_col./sum(norma_col);
Pcol = cumsum(pcol);
samp_col = @(k) nnz(rand>Pcol)+1;



%% perform RK method

maxiter = 1000;

rng(0)

x = zeros(size(A,2),1);

for iter = 1:maxiter
    
    r = samp_row(iter);
    %j = samp_col(iter);
    a = A(r,:);
    
    x = x - (a*x-b(r))/norma(r) *a';
    
end

disp(['Residual: ' num2str(norm(A*x-b))])
    



%% perform RSK method

maxiter = 10000;

[x,xdual] = deal( zeros(size(A,2),1) );

lambda = 0.5;
S = @(z) max(abs(z)-lambda,0).*sign(z); % Soft thresholding  

for iter = 1:maxiter
    
    r = samp_row(iter);
    a = A(r,:);
    xdual = xdual - (a*x-b(r))/norma(r)* a';
    x = S(xdual);
    
end

disp(['Residual: ' num2str(norm(A*x-b))])



%% perform REK method

maxiter = 1e7;

[x,xdual] = deal( zeros(size(A,2),1) );
z = b;

lambda = 0.5;
S = @(z) max(abs(z)-lambda,0).*sign(z); % Soft thresholding  

for iter = 1:maxiter
    
    r = samp_row(iter);
    s = samp_col(iter);
    a = A(r,:);
    c = A(:,s);
    
    z = z - (c.'*z)/norma_col(s)*conj(c);
    
    x = x - (a*x-b(r)+z(r))/norma(r)* a';
    %x = S(xdual);
    
    if mod(iter,1e4) == 0
        
        disp(['Distance to Moore-Penrose inverse solution: ' num2str(norm(x-pinv(A)*b))])
        
        %disp(['||ATz|| = ' num2str(norm(A.'*z))])
        
        proj_b = b - A*pinv(A)*b;
        disp(['||z_k-proj_b|| = ' num2str(norm(z-proj_b))])
        
        %disp(['Gradient of least-squares functional: ' num2str(norm(A.'*(A*x-b)))])
    end
    
end






%% perfrom some steps of ExRSK for comparison 

for iter = 1:4
    i = 2;
    j = 2;
    ai = A(i,:).';
    ajt = A(:,j);

    z = z - ajt.'*z/norm(ajt)^2 *ajt;
    xdual = xdual - (ai.'*x-b(i))/norm(ai)^2 *ai;
    x = S(xdual);
end

i = 1;
j = 2;
ai = A(i,:).';
ajt = A(:,j);

z = z - ajt.'*z/norm(ajt)^2 *ajt;
xdual = xdual - (ai.'*x-b(i))/norm(ai)^2 *ai;
x = S(xdual);

    