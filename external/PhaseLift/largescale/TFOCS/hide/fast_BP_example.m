% Make sure to add the solver path


randn('state',324324);

N       = 1e6;
sqrtN   = sqrt(N);
isqrtN  = 1/sqrtN;
M       = round(N/2);
rows    = randperm(N);
rows    = sort(rows(1:M));

FFT = linop_handles( [N,N], @(x) isqrtN*fft(x), @(x) sqrtN*real(ifft(x)) ,'R2C');
% linop_test(FFT)
A   = linop_compose( linop_subsample([M,N], rows ), FFT );
% linop_test(A)

% -- signal
K       = round(M/10);
x       = zeros(N,1);
T       = randperm(N);
T       = sort(T(1:K));
x(T)    = randn(K,1);
b       = A(x,1);

%% solve without continuation
mu      = 1;
x0      = 0;
z0      = [];
opts    = [];
opts.tol        = 1e-14;
opts.countOps   = 1;  % the "A" column is the # of FFT's
opts.printEvery = 10;
opts.errFcn     = {@(f,dual,primal) norm(primal - x )/norm(x), ...
                   @(f,dual,primal) norm(primal - x, Inf ) };
opts.restart    = 200;  % doesn't have too much effect
opts.maxIts     = 200;
               
% stop everything when we are at machine precision:
opts.stopFcn    = {@(f,dual,primal) (norm(primal-x,Inf) < 1e-15) };

[xk,out,optsOut]    = solver_sBP( A, b, mu, x0, z0, opts );

%% use continuation
mu          = 10;
opts.restart = 100; % doesn't have too much effect
contOpts    = [];
contOpts.innerMaxIts = 10;
contOpts.maxIts     = 4;
solver  = @(mu,x0,z0,opts )  solver_sBP( A, b, mu, x0, z0, opts );
tic
[xk,out,optsOut]   = continuation( solver, mu, x0, z0, opts, contOpts );
toc
%% plot
figure(1);
semilogy( out.err(:,1) );
hold all