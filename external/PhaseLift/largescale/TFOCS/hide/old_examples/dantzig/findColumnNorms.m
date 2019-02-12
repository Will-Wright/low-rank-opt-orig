function [d,Amatrix] = findColumnNorms(A,n,m)
% delta = findColumnNorms( A, N )
%   finds the norms of the columns of the matrix A
%   where A is given implicitly by a function handle
%
%   The memory used is O( max(M,N) )
%   as opposed to O( M*N ) if the full matrix were formed.
%
%   N is the dimension of the domain of A.
%
% delta = findColumnNorms( At, N, M )
%   returns the same result, but At should be a function handle
%   to the transpose matrix.  If M << N, this is faster.
%   (e.g. DCT with N = 16384, M=8192, then this calling sequence
%    takes about 5 seconds, whereas the first calling sequence takes 
%    49 seconds)
%
% [delta,Amatrix] = findColumnorms(...)
%   returns the explicit matrix for A (this IS O(M*N) memory)
%
% TFOCS version 1.0, code by Michael Grant (mcg@cvxr.com) and Stephen Becker (srbecker@caltech.edu)

if isa(A,'function_handle')
    % Now, have two modes: if m is provided, and m < n,
    % then assume that we really have been passed At,
    % and use this:
    if nargin > 2 && ~isempty(m) && isnumeric(m) && m < n && nargout < 2
        e = zeros(m,1);
        col = zeros(n,1);
        for k = 1:m
            e(k) = 1;
            col = col + A(e).^2;
            e(k) = 0;
        end
        d = sqrt( col )';
        
    else
    
        step = 1;  % step size, to minimze # of calls to A
        d = zeros(1,n);
        e = zeros(n,step);
        
        if nargout < 2
            if step > 1
                for k = 1:step:n
                    interval = [k:min(n,k+step-1)];
                    l = length(interval);
                    e(interval,1:l) = eye(l);
                    Ae = A(e);
                    d(interval) = sqrt(sum( Ae(:,1:l).^2 ) );
                    e(interval,1:l) = zeros(l);
                end
            else
                for k = 1:step:n
                    e(k) = 1; % the kth unit vector
                    d(k) = norm( A(e) );
                    e(k) = 0;
                end
            end
        else
            m = length(A(e));
            Amatrix = zeros(m,n);
            for k = 1:n
                e(k) = 1; % the kth unit vector
                Amatrix(:,k) = A(e);
                d(k) = norm( Amatrix(:,k) );
                e(k) = 0;
            end
        end
    end
else
    Amatrix = A;
    d = sqrt( sum(A.^2) );
end
