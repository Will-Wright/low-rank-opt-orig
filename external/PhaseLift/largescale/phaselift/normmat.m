function [AN,a] = normmat(A,flag,p);
% NORMMAT - computes column norms of A and normalizes the columns to 1
% Usage: [AN,a] = normmat(A,flag,p);
% Input: A    ... matrix
%        flag ... if flag = 0 norms of columns of A are computed but columns 
%                 are not normalized afterwards (output is vector a); 
%                 default: flag = 1.
%        p    ... p-norm is used, default: p = 2
% Output: AN  ... matrix with columns with l_2 norm 1
%          a  ... vector containing l_2 norms of columns of A
% 

% Author:     T. Strohmer, 2002 (strohmer@math.ucdavis.edu)
% Modified:
%
% Literature:  see my paper on Grassmannian frames
%
% Copyright:  (c) Thomas Strohmer, University of California, Davis, USA.
%             All standard disclaimers apply.
%


if nargin < 3
  p = 2;
end

if nargin < 2
  flag = 1;
end

[n,m]=size(A);
if flag == 0
  for k=1:m
    a(k)=norm(A(:,k),p);
  end
  AN=a;
else
  AN = zeros(size(A));
  for k=1:m
    a(k)=norm(A(:,k),p);
    AN(:,k)=A(:,k)/a(k);
  end
end
