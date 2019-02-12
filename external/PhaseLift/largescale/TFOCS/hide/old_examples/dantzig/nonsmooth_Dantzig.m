function [ val, x, approx ] = nonsmooth_Dantzig( varargin )

% Proximity operator for the l_1 norm ||x||_1
%   or the scaled l_1 norm || diag(s)*x ||_1
%
% One argument: 
%   val = L1_projector(x)    returns  val = sum(abs(x))
% Two arguments:
%   val = L1_projector(s,x)  returns  val = sum(abs(s.*x))
% Three and four arguments:
%   [val,x,approx] = L1_projector(y,g,L) (s=1)
%   [val,x,approx] = L1_projector(s,y,g,L)
%   Returns sum(abs(s.*x)) and the optimal point and value of
%       min_x g'*x + 0.5 * sum(L.*(x-y).^2) + sum(abs(s.*x))
%   "s" may be a vector of size(x) or a scalar.
%
% Warning: not yet tested with complex data

switch nargin,
    case 1,
        val = sum( abs( varargin{1} ) );
        return
    case 2,
        val = sum( abs( varargin{1} .* varargin{2} ) );
        return
    case 3,
        s = 1;
    case 4,
        s = varargin{1};
end
[ y, g, L ] = deal( varargin{nargin-2:end} );
if ~isreal(y) || ~isreal(g) || ~isreal(L) || ~isreal(s)
    % --- Complex numbers case ---
    % This is a rough first-pass
    x = y - g./L;
    scale = max( 0, 1 - s./( L.*abs(y) ) );  % unlike prox_L2,
                                             % this is a vector, not a scalar
	x = x .* scale;
    val = sum(abs(s.*x));
    if nargout > 2
        approx = g'*(x-y) + 0.5*sum(L.*(x-y).^2) + sum(abs(s.*x));
    end
else
    % --- Real numbers case --
    
    % x = shrink( y - g ./ L, abs(s) ./ L )
    x = L .* y - g;
    x = sign( x ) .* max( abs(x) - abs(s), 0 );
    x = x ./ ( L + realmin * ( x == 0 ) );
    val = sum(abs(s.*x));
    if nargout > 2,
        % For the line below, we want g'*x, not g'*(x-y), no?  At least that's
        % what the help text at the beginning of the file suggests.
        
        % g'*(x-y) + 0.5*sum(L.*(x-y).^2) + sum(abs(s.*x));
        
        % Same as above line, but more stable?
        approx = -sum(L.*(x-y).*(x+y))/2;
    end
end