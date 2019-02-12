function Ay = linear_Dantzig( A, At, b, D, lambda, adjoint )
% Ay = linear_Dantzig( A, At, b, D, lambda, adjoint )
%   represents a linear operator A.
%   Setting ADJOINT = 'adjoint' will compute A^*
%       otherwise it will compute A

% last edit: April 14, 2010

% a little bit of extra overhead, but simpler for the user:
if ~isa(At,'function_handle')
    if isempty(At)
        if isa(A,'function_handle')
            error('When not providing At, A must be a matrix');
        end
        At = @(x)A'*x;
    else
        At = @(x)At*x;
    end
end
if ~isa(A,'function_handle')
    A = @(x)A*x;
end
if ~isa(D,'function_handle')
    if ~isvector(D)
        D = diag(D);
    end
    D = @(x) D.*x;
end

error(nargchk(5,6,nargin));
MODE = 1; % default.  Use all vectors to make everything linear
% MODE = 2; % new.  Use cells to contain both linear and offset part
if nargin == 6 && strcmpi( adjoint, 'adjoint' )
    % (transpose) (lambda is really dual_x )
    if MODE == 1
        Ay = A( lambda(1:end-1) ) + b * lambda(end);
    else
        Ay = A( lambda ) - b;
    end
    Ay = D( At( Ay ) );
else
    % (forward)
    Ay = A( D( lambda ) );
    if MODE == 1
        Ay = [ At( Ay ) ; b' * Ay ];
    else
        Ay = { At(Ay), b'*Ay };
    end
end
