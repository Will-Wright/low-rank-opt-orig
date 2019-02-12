function [evecs,evals] = rightmostEigenproblem(A,y,opts,v)
   if nargin == 4 && ~isempty(v)
      eigsopts.v = v;
   end
   if nargin < 3 || isempty(opts)
      opts = struct('isreal',false,'eigsk',1,'eigsp',3,'eigstol',double(eps('single')));
   end
   if ~isfield(opts,'isreal') || isempty(opts.isreal)
      opts.isreal = false;
   end
   if ~isfield(opts,'eigsk') || isempty(opts.eigsk)
      opts.eigsk = 1;
   end
   if ~isfield(opts,'eigsp') || isempty(opts.eigsp)
      opts.eigsp = 3;
   end
   if ~isfield(opts,'eigstol') || isempty(opts.eigstol)
      opts.eigstol = eps;
   end
   if opts.isreal
      vec   = @(x)real(x(:));
      sigma = 'LA';
   else
      vec   = @(x)x(:);
      sigma = 'LR';
   end
   Acty = A.adjoint(y);
   Afun = @(x)vec(Acty(x));
   n    = size(A,2);
   k    = opts.eigsk;
   eigsopts.isreal     = opts.isreal;
   eigsopts.issym      = true;
   eigsopts.p          = opts.eigsp;
   eigsopts.tol        = opts.eigstol;
   [evecs,evals]       = eigs(Afun,n,k,sigma,eigsopts);
   [evals,permutation] = sort(real(diag(evals)),'descend');
   evecs               = evecs(:,permutation);
end
