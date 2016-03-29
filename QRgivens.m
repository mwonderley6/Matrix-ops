
function [x,rhoLS] = QRgivens(A,b);
%
%  Input:   A  -  an m by n matrix with m .ge. n
%
%           b  -  an m vector
%
%  Output:  x  -  The solution to min norm(b - Ax)
%
%           rhoLS = min norm(b - Ax) 
%
%  Factors A = QR ,  where A is m by n
%  Q'Q = I_m, R is m by n upper triangular
%  by Givens method and reduces b at the 
%  same time.
%
%  Then the solution  to Rx = b1 is obtained
%
%
%  D.C. Sorensen 
%  30 Oct 00
%-----------------------------------------------------------------

   [m,n] = size(A);
   tau = zeros(n,1);
   if (m == n), x = A\b; rhoLS = 0; return; end
%
%  Reduce the initial n+1 by n+1 submatrix of [A b]
%  to upper triangular form
%
   R = [A(1:n+1,:) b(1:n+1)] ;
   for j = 2:n+1,
      for i = 1:j-1,
         [c,s,r] = Givens(R(i,i),R(j,i));
         R(i,i) = r; R(j,i) = 0;
         t =             c*R(i, i+1:n+1) + s*R(j,i+1:n+1);
         R(j,i+1:n+1) = -s*R(i, i+1:n+1) + c*R(j,i+1:n+1);
         R(i, i+1:n+1) = t;
      end
   end
%
%  Process the remaining m - (n+1) rows of [A b]
%
   for k = n+2:m,
      a = [A(k,:) b(k)];
      for i = 1:n+1,
         [c,s,r] = Givens(R(i,i),a(i));
         R(i,i) = r; a(i) = 0;
         t =           c*R(i, i+1:n+1) + s*a(i+1:n+1);
         a(i+1:n+1) = -s*R(i, i+1:n+1) + c*a(i+1:n+1);
         R(i, i+1:n+1) = t;
      end
   end
%
%  Solve the upper triangular system  Rx = b1 in place
%  b1 resides in R(1:n,n+1);
%
   x = R(1:n,n+1);
   for j = n:-1:2
       x(j) = x(j)/R(j,j);
       x(1:j-1) = x(1:j-1) - R(1:j-1,j)*x(j);
   end 
   x(1) = x(1)/R(1,1);
   rhoLS = abs(R(n+1,n+1));
