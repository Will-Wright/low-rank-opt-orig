function [U,S,V,flag] = svds(varargin)
% SVDS for packSVD class
%   See the documentation for the builtin SVDS
%
%   This calls builtin EIGS
%
%   See also SVDS, EIGS


A = varargin{1};
[m,n] = size(A);

p = min(m,n);
q = max(m,n);

% B's positive eigenvalues are the singular values of A
% "top" of B's eigenvectors correspond to left singular values of A
% "bottom" of B's eigenvectors correspond to right singular vectors of A

% B = [sparse(m,m) A; A' sparse(n,n)];
B = @(x) [ A*x(m+1:end); A'*x(1:m) ];   % SRB
sz = m + n;

% Overload the normest operator:
Af = @(x) A'*(A*x);
normest = @(X) eigs(Af, n, 1);  % SRZB

if nargin < 2
   k = min(p,6);
else
   k = varargin{2};
end
if nnz(A) == 0
   if nargout <= 1
      U = zeros(k,1);
   else
      U = eye(m,k);
      S = zeros(k,k);
      V = eye(n,k);
      flag = 0;
   end
   return
end
if nargin < 3
   bk = min(p,k);
   if isreal(A)
       bsigma = 'LA';
   else
       bsigma = 'LR';
   end
else
   sigma = varargin{3};
   if sigma == 0 % compute a few extra eigenvalues to be safe
      bk = 2 * min(p,k);
   else
      bk = k;
   end
   if strcmp(sigma,'L')
       if isreal(A)
           bsigma = 'LA';
       else
           bsigma = 'LR';
       end
   elseif isa(sigma,'double')
      bsigma = sigma;
      if ~isreal(bsigma)
          error('MATLAB:svds:ComplexSigma', 'Sigma must be real.');
      end
   else
      error('MATLAB:svds:InvalidArg3',...
            'Third argument must be a scalar or the string ''L''.')
   end   
end
if nargin < 4
   % norm(B*W-W*D,1) / norm(B,1) <= tol / sqrt(2)
   % => norm(A*V-U*S,1) / norm(A,1) <= tol
   boptions.tol = 1e-10 / sqrt(2);
   boptions.disp = 0;
else
   options = varargin{4};
   if isstruct(options)
      if isfield(options,'tol')
         boptions.tol = options.tol / sqrt(2);
      else
         boptions.tol = 1e-10 / sqrt(2);
      end
      if isfield(options,'maxit')
         boptions.maxit = options.maxit;
      end
      if isfield(options,'disp')
         boptions.disp = options.disp;
      else
         boptions.disp = 0;
      end
   else
      error('MATLAB:svds:Arg4NotOptionsStruct',...
            'Fourth argument must be a structure of options.')
   end   
end

% [W,D,bflag] = eigs(B,bk,bsigma,boptions);
boptions.isreal = isreal(A); % SRB
boptions.issym = true; % SRB
[W,D,bflag] = eigs(B,sz,bk,bsigma,boptions);  % SRB

if ~isreal(D) || (isreal(A) && ~isreal(W))
   error('MATLAB:svds:ComplexValuesReturned',...
         ['eigs([0 A; A'' 0]) returned complex values -' ...
         ' singular values of A cannot be computed in this way.'])
end

% permute D and W so diagonal entries of D are sorted by proximity to sigma
d = diag(D);
if strcmp(bsigma,'LA') || strcmp(bsigma,'LR')
   nA = max(d);
   [dum,ind] = sort(nA-d);
else
   nA = normest(A);
   [dum,ind] = sort(abs(d-bsigma));
end
d = d(ind);
W = W(:,ind);

% Tolerance to determine the "small" singular values of A.
% If eigs did not converge, give extra leeway.
if bflag 
   dtol = q * nA * sqrt(eps);
   uvtol = m * sqrt(sqrt(eps));
else
   dtol = q * nA * eps;
   uvtol = m * sqrt(eps);
end

% Which (left singular) vectors are already orthogonal, with norm 1/sqrt(2)?
UU = W(1:m,:)' * W(1:m,:);
dUU = diag(UU);
VV = W(m+(1:n),:)' * W(m+(1:n),:);
dVV = diag(VV);
indpos = find((d > dtol) & (abs(dUU-0.5) <= uvtol) & (abs(dVV-0.5) <= uvtol));
indpos = indpos(1:min(end,k));
npos = length(indpos);
U = sqrt(2) * W(1:m,indpos);
s = d(indpos);
V = sqrt(2) * W(m+(1:n),indpos);

% There may be 2*(p-rank(A)) zero eigenvalues of B corresponding
% to the rank deficiency of A and up to q-p zero eigenvalues
% of B corresponding to the difference between m and n.

if npos < k
   indzero = find(abs(d) <= dtol);
   QWU = orth(W(1:m,indzero));
   QWV = orth(W(m+(1:n),indzero));
   nzero = min([size(QWU,2), size(QWV,2), k-npos]);
   U = [U QWU(:,1:nzero)];
   s = [s; abs(d(indzero(1:nzero)))];
   V = [V QWV(:,1:nzero)];
end

% sort the singular values in descending order (as in svd)
[s,ind] = sort(s);
s = s(end:-1:1);
if nargout <= 1
   U = s;
else
   U = U(:,ind(end:-1:1));
   S = diag(s);
   V = V(:,ind(end:-1:1));
%    flag = norm(A*V-U*S,1) > sqrt(2) * boptions.tol * norm(A,1);
   flag = norm(A*V-U*S,1) > sqrt(2) * boptions.tol * norm(A,'fro');  % SRB
end


