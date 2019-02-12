% Ewout's test, modified by Stephen to test a more complicated case,
%  and with complex numbers.

% Solve:  min_X  mu/2||X-Y||_F^2     s.t.  X >= 0


% ----------------------------------------------
smooth = smooth_quad;
proj   = @proj_psd;

N = 11;

xInit = randn(N);
xInit = xInit + 1i*randn(N);
xInit = xInit * xInit';

Y       = randn(N);     % not PSD
Y       = Y + 1i*randn(N);
Y       = Y + Y';       % but Y is at least Hermitian

vec     = @(X) X(:);
mat     = @(x) reshape(x,N,N);
A       = linop_handles({[N,N],[N^2,1]}, vec, mat,'c2c');       % identity

% -- starting point --
% Does this have to be Hermitian?  Apparently not
% x0      = xInit;
% x0      = Y;
% x0      = randn(N);
x0      = randn(N) + 1i*randn(N);

% Compute the exact answer, which we use as a reference
[V,D]   = eig(Y);
d       = diag(D);
pos     = find( d > 0 );
X_answer = V(:,pos)*diag(d(pos))*V(:,pos)';

opts = [];
% opts.alg        = 'AT';
% otps.restart    = 100;
opts.alg        = 'GRA';    % works very well
opts.maxIts     = 1000;
opts.tol        = 1e-13;
opts.printEvery = 10;
opts.maxmin     = 1; % Minimize
% opts.L0         = 1000;
% opts.beta       = 1;

opts.errFcn     = @(f,X) norm(X-X_answer,'fro');

% also fails if xInit is zero.  Bad!
[x,out,opts] = tfocs(smooth,{A,-vec(Y)},proj,x0, opts);
% ----------------------------------------------


% -- Check our answer --
% x should be the PSD projection of Y, i.e. just the non-negative
% eigenvalues.  We can check that:
disp(sort([eig(x), eig(Y), real(eig(X_answer)) ],'descend'))