function out=ComputeErrors( xTrue, fTrue, x, f )
% ComputeErrors( xTrue, fTrue, x, f )
%   computes several errors relating the variable x to xTrue
%   and the objective function value f to fTrue.
%
%   There is no output; rather, a global variable CURRENT_ERROR is 
%   updated.  This is intended to be used with multiplyA, and record
%   the error at any given # of calls to A and At.
%   If an output is requested, the output is set to 0.
%
%   The variable CURRENT_ERROR may either be a scalar or a row-vector
%   (modify the code in this file to suit your needs).
%
%   TFOCS version 1.0, code by Michael Grant (mcg@cvxr.com) and Stephen Becker (srbecker@caltech.edu)
%
%   See also multiplyA

global CURRENT_ERROR

% Update this, so that other functions don't need to know the size
if nargin <= 2
    e = Inf( 1, 4 );
else

    e(1)    = norm(x-xTrue)/norm(xTrue);
    e(end+1)= norm(x-xTrue,'inf');
    
    T1 = find(xTrue);
    T2 = find(x);
    e(end+1)= length(T1)+length(T2)-2*length(intersect(T1,T2) );
    
%     e(end+1)= abs(f-fTrue)/abs(fTrue);
    e(end+1)= abs( norm(x,1) -fTrue)/abs(fTrue);


end


CURRENT_ERROR = e;
if nargout > 0
    out = e(2);
end
