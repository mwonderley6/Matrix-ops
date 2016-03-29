function [QR,tau] = QRfac(A);
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

   [m,n] = size(A);
   tau = zeros(n,1);

   for k = 1:n,
%
%    Compute the Householder vector v
%
%    I - tau(k)*v*v' is the transformation 
%
     v = A(k:m,k);
     rho = sign(v(1))*norm(v);
     if (abs(rho) > 0),
        v(1) = v(1) + rho;
        tau(k) = v(1)/rho;
        v = v/v(1);

        for j = k+1:n
%
%           Apply the Householder transformation to the j-th column of A
%
            alpha = tau(k)*(v'*A(k:m,j));
            A(k:m,j) = A(k:m,j) - v*alpha;

        end
        A(k:m,k) = v;
        A(k,k) = -rho;

     end
   end
   QR = A;
     


   

