function varargout = size( x, dim )

% SIZE   TFOCS-friendly size operator.

m = x.sz(1);
n = x.sz(2);

if nargin > 1
    if dim == 1
        varargout{1} = m;
    elseif dim == 2
        varargout{1} = n;
    else
        error('bad value for dimension');
    end
    return;
end

if nargout < 2
    varargout{1} = [ m, n ];
else
    varargout{1} = m;
    varargout{2} = n;
end

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
