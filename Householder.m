function[Q,R] = Householder(A)
[m,n] = size(A);
U = A(:,1);
R = A;
fprintf('U = %s\n',num2str(U));
fprintf('R = %s\n',num2str(R));

%fprintf('x=%s',num2str(R));

for i = 1:m
    fprintf('I = %s\n',num2str(i));
    U = A((((i-1)*m)+i):(i*m));
fprintf('Ui = %s\n',num2str(U));
    normU = norm(U);
fprintf('normU = %s\n',num2str(normU));
    U(1) = U(1) + normU;
    fprintf('x = %s\n',num2str(U));
    Q = eye(m,n)-((-2*(U*U'))/(norm(U)^2));
    R = Q*R;
    fprintf('x = %s\n',num2str(Q));
end

end