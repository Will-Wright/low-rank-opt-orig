function ans = nnz( x )
ans = sum( cellfun( @nnz, x.value_ ) );

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
