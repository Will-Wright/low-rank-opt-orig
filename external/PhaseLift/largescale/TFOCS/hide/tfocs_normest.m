function nrm = tfocs_normest( A, squared, transpose_string )
% norm = tfocs_normest(A)
%   estimates the norm of A
%
% norm = tfocs_normest(A,2)
%   estimates the norm of A*A'
%
% norm = tfocs_normest(A,2,'transpose')
%   estimates the norm of A'*A
%
% "A" is either a matrix or a TFOCS operator

error(nargchk(1,3,nargin));
if nargin < 2 || isempty(squared)
    squared = false;
else
    squared = true;
    transposed = false;
    if nargin > 3 && (...
        (ischar(transpose_string) && ~isempty(strfind(lower(transpose_string),'tran')) ) ...
        || (~ischar(transpose_string) && transpose_string )  )
        transposed = true;
    end
end

tol     = 1e-5;
maxiter = 100;

if isnumeric(A)
    if squared && transposed
        A = A'*A;
    elseif squared
        A = A*A';
    end
    if numel(A) < 1000^2
        nrm = norm(A);
    else
        nrm = my_normest(A,tol,maxiter);
    end
elseif isa(A,'function_handle')
    
    sz = A([],0);
    if iscell(sz)
         n = sz{1};
         m = sz{2};
    else
        m = sz(1);
        n = sz(2);
    end
    
    Size = n;
    if squared && transposed    % A'*A
        Af = @(x) A(A(x,1),2);
        At = Af;
    elseif squared              % A*A'
        Af = @(x) A(A(x,2),1);
        At = Af;
        Size = m;
    else                        % A
        Af = @(x) A(x,1);
        At = @(x) A(x,2);
    end
    
    if max(Size) > 1e4
        tol = 1e-4;
        maxiter = 50;
    end
    if max(Size) > 1e5
        tol = 1e-3;
        maxiter = 25;
    end
    nrm = my_normest(Af,At,Size,tol,maxiter);
    
else
    error('A must be a matrix or function handle');
end
        
% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.