function g = mrdivide(a,b)
% handles W/alpha, where W is a function, alpha a scalar

if isscalar(b) && isa(a,'function_handle')
    g = mtimes(a,1/b);
    return;
end
error('MATLAB:UndefinedFunction','Undefined function or method %s for input arguments of type %s and %s',mfilename,class(a),class(b) );

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
