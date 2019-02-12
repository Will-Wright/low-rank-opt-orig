function g = mtimes(a,b)

if isscalar(a) && isa(b,'function_handle')
    if isProx( func2str(b) )
        g = prox_scale(b,a);
        return
    end
    if isLinear( func2str(b) )
        g = linop_compose( linop_scale(a), b );
        return
    end
elseif isscalar(b) && isa(a,'function_handle')
    if isProx( func2str(a) )
        g = prox_scale(a,b);
        return
    end
    if isLinear( func2str(a) )
        g = linop_compose( linop_scale(b), a );
        return
    end
end

error('MATLAB:UndefinedFunction','Undefined function or method %s for input arguments of type %s and %s',mfilename,class(a),class(b) );


function good = isProx( f )
% We could explicitly make a list, or for now, we just
% assume that all the proximity functions are named
% like "prox_****"
good = false;
if strfind(f,'prox_') == 1
    good = true;
end

function good = isLinear( f )
% We could explicitly make a list, or for now, we just
% assume that all the proximity functions are named
% like "linop_****"
good = false;
s = strfind(f,'linop_');
if s == 1
good = true;
elseif ~isempty(s) && s>3 && all(f(1:2) == '@(') && f(s(1)-1) == ')'
    % this catches things like 
    % f ='@(x,mode)linop_matrix_r2r(sz,A,x,mode)'
    good = true;
end

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
