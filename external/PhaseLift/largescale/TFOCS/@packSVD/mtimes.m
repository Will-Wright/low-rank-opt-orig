function z = mtimes( x, y )

% MTIMES   Multiplication. Multiplication by scalars only.

% z = times( x, y );
if isscalar(y)
    % we scale the matrix, rather than do a matrix multiply
    z = x;
    if ~isempty(x.X)
        z.X = z.X*y;
    end
    if ~isempty(x.s)
        z.s = z.s*y;
    end
    return;
end
if ~isa(x,'packSVD') && isscalar(x)
    z = y;
    if ~isempty(y.X)
        z.X = y.X*x;
    end
    if ~isempty(y.s)
        z.s = y.s*x;
    end
    return;
end

% z = zeros( size(x,1), size(y,2) );
z1 = 0;
z2 = 0;
if ~isempty(x.X)
    z1 = x.X*y;
end
if ~isempty(x.s)
    z2 = x.U*( diag(x.s) * ( x.V'*y) );
end
z = z1 + z2;



% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
