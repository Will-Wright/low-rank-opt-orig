function [e,cnt] = my_normest(S,St,n,tol, maxiter)
%NORMEST Estimate the matrix 2-norm.
%   NORMEST(S) is an estimate of the 2-norm of the matrix S.
%   NORMEST(S,tol) uses relative error tol instead of 1.e-6.
%   [nrm,cnt] = NORMEST(..) also gives the number of iterations used.
%
%   NORMEST(S,St,n,[tol]) uses the funcions S() and St() to compute
%   the forward and transpose multiplication (where S is m x n)
%   This modification due to Stephen Becker, 10/01/09
%   If St is the empty matrix, then we assume S = S'
%   Note: this method is just the power method.
%
%   This function is intended primarily for sparse matrices,
%   although it works correctly and may be useful for large, full
%   matrices as well.  Use NORMEST when your problem is large
%   enough that NORM takes too long to compute and an approximate
%   norm is acceptable.
%
% Version modified by Stephen Becker to allow for function handles
% as part of
% TFOCS version 1.0, code by Michael Grant (mcg@cvxr.com) and Stephen Becker (srbecker@caltech.edu)
%
%   See also normest

%   Copyright 1984-2006 The MathWorks, Inc. 
%   $Revision: 5.14.4.3 $  $Date: 2006/10/02 16:32:45 $

if nargin < 2, tol = 1.e-6; end
if nargin < 5, maxiter = 20; end
IMPLICIT = false;
if isa(S,'function_handle')
    if isempty(St)
        St = S;  % we assume the matrix is symmetric;
    elseif ~isa(St,'function_handle')
        error('normest: must provide transpose function');
    end
    if nargin < 3
        error('normest: must provide width of matrix');
    end
    
    if nargin < 4, tol = 1.e-6; end
    IMPLICIT = true;
else
    if nargin < 3 && isnumeric(St), tol = St; end
end

if ~IMPLICIT
    x = sum(abs(S),1)';
else
    x = ones(n,1);
end

cnt = 0;
e = norm(x);
if e == 0, return, end
x = x/e;
e0 = 0;
while abs(e-e0) > tol*e && cnt < maxiter
   e0 = e;
   if ~IMPLICIT
       Sx = S*x;
   else
       Sx = S(x);
   end
   if nnz(Sx) == 0
      Sx = rand(size(Sx));
   end
   e = norm(Sx);
   if ~IMPLICIT
       x = S'*Sx;
   else
       x = St(Sx);
   end
   x = x/norm(x);
   cnt = cnt+1;
end
