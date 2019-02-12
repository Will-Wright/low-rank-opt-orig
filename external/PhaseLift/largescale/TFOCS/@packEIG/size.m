function v = size( x )

% SIZE   TFOCS-friendly size operator.

m = x.sz(1);
n = x.sz(2);


v = [ m, n ];

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
