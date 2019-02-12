function [ v, U ] = proj_tsvd( ApplyZhandle, n, RNK )

% [ v, U ] = proj_tsvd( ApplyZhandle, n , RNK)
% Implements the truncated SVD projection keeping at most RNK singular values
% Assumes that the matrix is hermitian
% ApplyZhandle is a handle implicitly computing the gradient according to 
% AugNegLoglikelihood

v = 0;

eigsopts.issym = true;
eigsopts.isreal = false;

% MPF Jul 2, 2015: Change 'LM' to 'LR'.
[V,D] = eigs(ApplyZhandle,n,RNK,'LR',eigsopts);

D = max(diag(D),0);
tt = D > 0;
tt = 1:length(D);
U  = bsxfun(@times,V(:,tt),sqrt(D(tt,:))');

if min(D) < -8*eps*max(D),
    v = Inf;
end


