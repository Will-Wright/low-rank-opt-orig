function op = linop_scale( scale )

%LINOP_DOT  Scaling linear operator.
%    OP = LINOP_SCALE( scale ) returns a handle to a TFOCS linear operator 
%    whose forward and adjoint operators are OP(X) = scale * X.

if ~isnumeric( scale ) && numel( scale ) ~= 1,
    error( 'Argument must be a scalar.' );
elseif ~isreal( scale ),
    error( 'Argument must be real.' );
end
if scale == 1,
    op = @linop_identity;
else
    op = @(x,mode)linop_scale_impl( scale, x, mode );
end

function y = linop_scale_impl( scale, y, mode )
if mode == 0, 
    y = { [], [] };
else
    y = scale * y;
end

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
