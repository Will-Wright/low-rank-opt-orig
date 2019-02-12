function op = smooth_logdet()
% Function to represent -log( det(x) )
%   (Note the minus sign)
% x must be Hermitian and positive definite.

% SRB: have not yet tested this.
% SRB: I think we CAN compute the proximity operator to logdet.
%       Will implement this in prox_logdet
op = @smooth_logdet_impl;

function [ v, g ] = smooth_logdet_impl( x )
if size(x,1) ~= size(x,2)
    error('smooth_logdet: input must be a square matrix');
end
if nnz(x) == 0
    % Last step hopefully
    v = - Inf;
    if nargout > 1, g = x; end
    return;
end
v = -log(det(x));
% Note: we may wish to calculate ourself, e.g.:
%R = chol(x);
%v = -sum(log(diag(R)));

if nargout > 1,
    % Is there a minus sign here?  Need to check
    g = -inv(x);
    % it would be nice to be able to return a function handle to g,
    % and then calculate g(y) as x\y
end

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
