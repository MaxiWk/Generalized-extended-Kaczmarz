
A = [1,2,3;
     4,5,6];
xhat = [2,0,1]';
b = A*xhat;

epsilon = 0.1;
lambda = 0.1;
L_gstar = max(1/epsilon,1);


zdual = b;
z = T_repsilon(zdual,epsilon);
[x, xdual] = deal( zeros(size(xhat)) );
     
j = 3;
i = 1;
zdual = zdual - A(:,j)'*z/(L_gstar*norm(A(:,j))^2)* A(:,j);
z = T_repsilon(z,epsilon);

xdual = xdual - (A(i,:)*x - b(i) + zdual(i))/norm(A(i,:))^2 * A(i,:)';
x = max(abs(xdual)-lambda,0).*sign(xdual); 




j = 2;
i = 2;
zdual = zdual - A(:,j)'*z/(L_gstar*norm(A(:,j))^2)* A(:,j);
z = T_repsilon(z,epsilon);

xdual = xdual - (A(i,:)*x - b(i) + zdual(i))/norm(A(i,:))^2 * A(i,:)';
x = max(abs(xdual)-lambda,0).*sign(xdual); 




j = 2;
i = 2;
zdual = zdual - A(:,j)'*z/(L_gstar*norm(A(:,j))^2)* A(:,j);
z = T_repsilon(z,epsilon);

xdual = xdual - (A(i,:)*x - b(i) + zdual(i))/norm(A(i,:))^2 * A(i,:)';
x = max(abs(xdual)-lambda,0).*sign(xdual); 