function [ varargout ] = eigs( x, k )
% D = EIGS( X )
% [V,D] = EIGS( X )
% ... = EIGS( X, k) returns k largest eigenvalues
%
% X is a packSVD object with fields U, s and V
% Only the V and s fields are used.
%   
% See also eigs

if nargin < 2, k = 6; end
% EIGS   EIGS.

if x.orth && issparse( x.X ) && ~nnz( x.X ),

    s = abs( x.s );
    if nargout > 1,
	    U = x.U; 
    	V = conj(x.V);
    	if ( ~isreal( s ) || any( s < 0 ) ),
    		if x.sz(1) < x.sz(2),
	        	U = bsxfun( @times, U, x.s ./ s );
	        else
		        V = bsxfun( @times, V, x.s ./ s );
			end
		end	
    end
    
elseif nargout == 1,

    % we can as for 'econ' here, can't we??
	%s = eigs(full(double(x)));
	s = eigs(full(double(x)),'econ');
	
else	

    [U,s,V] = eigs(full(double(x)),'econ');
    s  = diag(s)';
    tt = s == 0;
    if any(tt),
    	s(:,tt) = [];
    	U(:,tt) = [];
    	V(:,tt) = [];
	end
	
end

if nargout >= 2
    varargout{1} = U;
    varargout{2} = diag(s);
    if nargout >= 3
        varargout{3} = V;
    end
else
    varargout{1} = s;
end


% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
