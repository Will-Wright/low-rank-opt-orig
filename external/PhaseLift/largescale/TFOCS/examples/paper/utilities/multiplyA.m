function out = multiplyA( A, At, type, x )
% out = multiplyA( A, At, type, x )
%   should be fashioned into "forward" and "transpose/adjoint" calls as follows:
%
%   Aforward   = @(x) multiplyA( A, At, 'forward', x )
%   Atranspose = @(x) multiplyA( A, At, 'transpose', x )
%
% where "A" is either a matrix (in which case "At" is ignored)
%   or a function handle (in which case "At" is the transpose function).
%
% The "type" can either be 'forward' (default) or 'transpose' or 'adjoint'
% for the transpose operation.
% Note: no distinction between 'transpose' and 'adjoint' is made.
%   It always assumes that the adjoint is wanted.  For real numbers,
%   this doesn't matter. 
%   (If "A" is a matrix, then out = A'*x )
%
% The purpose of this function is to update a persistent variable ERROR
% that contains information on the current error at a given iteration
% (where the iterations are measured in number of calls to A and At).
% The first column of this variable is the # of calls to A and At,
% and the second (and greater) column(s) is/(are) the error at that step.
%
% To access ERROR, call this function without any arguments, like:
%   
%   ERROR = multiplyA();
%
% Note: ERROR may be a row-vector, not just a scalar.
%
% This will also reset the value of the variable.
%
% To get the value of ERROR without resetting, call as
%   ERROR = multiplyA('noreset');
%
% TFOCS version 1.0, code by Michael Grant (mcg@cvxr.com) and Stephen Becker (srbecker@caltech.edu)
%
%   See also ComputeErrors

persistent ERROR
global CURRENT_ERROR

if nargin == 0
    out = ERROR;
    ERROR = [];
    % Also, this is important -- otherwise CURRENT_ERROR isn't reset
    ComputeErrors();
    return;
end
if nargin == 1 && ischar(A) && strcmpi(A,'noreset')
    out = ERROR;
    return;
end

if ~isempty(type) && ( ~isempty(strfind( lower(type), 'trans')) || ~isempty(strfind(lower(type),'adj')) )
    % Compute adjoint/transpose
    if isa( At, 'function_handle' )
        out = At(x);
    else
        out = A'*x;
    end
else
    if isa( A, 'function_handle' )
        out = A(x);
    else
        out = A*x;
    end
end
% update counter
if isempty(ERROR)
    ERROR = [1, CURRENT_ERROR];
else
    iter = ERROR(end,1);
    ERROR = [ ERROR; [iter+1, CURRENT_ERROR] ];
end
