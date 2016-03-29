
function [QR,tau] = QRfac2(A);
%
%  Input:   A  -  an m by n matrix with m .ge. n
%
%  Output:  QR -  The Householder QR factorization
%                 in product form.  QR can overwrite A.
%
%           tau - A vector of length n containing 
%                 "pivot" information
%
%  Factors A = QR ,  where A is m by n
%  Q'Q = I_m, R is m by n upper triangular
%
%  The leading n by n block of R is an 
%  n by n upper triangular matrix R_n
%
%  In practic A is overwritten by the 
%  Q in product form and R_n  in the 
%  upper triangle
%
%  D.C. Sorensen 
%  3 Oct 00
%-----------------------------------------------------------------
%
   [m,n] = size(A);
   tau = zeros(n,1);

   for k = 1:n,
%
%    Compute the Householder vector v
%
%    I - tau(k)*v*v' is the transformation  
%                    v = [1 ; A(k+1:m,k)];
%
     rho = sign(A(k,k))*norm(A(k:m,k));
     if (abs(rho) > 0),
        A(k,k) = A(k,k) + rho;
        tau(k) = A(k,k)/rho;
        A(k:m,k) = A(k:m,k)/A(k,k);
   
        for j = k+1:n
%
%           Apply the Householder transformation to the j-th column of A
%
            alpha = tau(k)*(A(k:m,k)'*A(k:m,j));
            A(k:m,j) = A(k:m,j) - A(k:m,k)*alpha;

        end
        A(k,k) = -rho;
   
     end
   end
   QR = A;
     


   

