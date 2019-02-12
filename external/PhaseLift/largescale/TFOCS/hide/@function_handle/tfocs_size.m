function g = size(a)

if isLinear( func2str(a) )
    if nargout == 0
        disp('TFOCS_operator of size');
    end
    sz = a([],0);
    if iscell(sz ) && nargout == 0
        fprintf(' [%d , %d ] x [%d , %d ]\n',...
            sz{2}(1), sz{2}(2), sz{1}(1), sz{1}(2) );
        return;
    end
    g = sz;
else
    g=builtin('size',a);
end


%error('MATLAB:UndefinedFunction','Undefined function or method %s for input arguments of type %s',mfilename,class(a) );


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
