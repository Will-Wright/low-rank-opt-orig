%{
    Tests basis pursuit

    min_x ||X||_1
s.t.
    A(X) == b

The solvers solve a regularized version, using
||x||_1 + mu/2*||x-x_0||_2^2

%}

% Before running this, please add the TFOCS base directory to your path

% Try to load the problem from disk
fileName = 'basispursuit_problem1_smoothed_noiseless';
randn('state',34324);

% We don't want to store "A"
rand('state',34324);
N = 1024;
M = round(N/2);
K = round(M/5);
A = randn(M,N);
if exist([fileName,'.mat'],'file')
    load(fileName);
    fprintf('Loaded problem from %s\n', fileName );
else
    % Generate a new problem
    x = zeros(N,1);
    T = randsample(N,K);
    x(T) = randn(K,1);

    b = A*x;
    EPS = 0;
    b_original = b;
    x_original = x;
    
    mu = .01*norm(x,Inf);
    x0 = zeros(N,1);
    % Note: with equality constraints, this is an LP, so for mu small
    % enough (but > 0), we have exact relaxation.

    % get reference via CVX
    tic
    cvx_begin
        cvx_precision best
        variable xcvx(N,1)
        minimize norm(xcvx,1) + mu/2*sum_square(xcvx-x0)
        subject to
            A*xcvx == b
    cvx_end
    time_IPM = toc;
    x_ref = xcvx;      
    obj_ref = norm(x_ref,1) + mu/2*sum_square(x_ref-x0);
    
    save(fileName,'x_ref','b','x_original','mu',...
        'EPS','b_original','obj_ref','x0','time_IPM');
    fprintf('Saved data to file %s\n', fileName);
    
end


[M,N]           = size(A);
K               = nnz(x_original);
norm_x_ref      = norm(x_ref);
norm_x_orig     = norm(x_original);
er_ref          = @(x) norm(x-x_ref)/norm_x_ref;
er_signal       = @(x) norm(x-x_original)/norm_x_orig;
resid           = @(x) norm(A*x-b)/norm(b);  % change if b is noisy

fprintf('\tA is %d x %d, original signal has %d nonzeros\n', M, N, K );
fprintf('\tl1-norm solution and original signal differ by %.2e (mu = %.2e)\n', ...
    norm(x_ref - x_original)/norm(x_original),mu );

%% Call the TFOCS solver
% er              = er_ref;  % error with reference solution (from IPM)
er              = er_signal; % error from original signal
opts = [];
opts.errFcn     = { @(f,dual,primal) er(primal), ...
                    @(f,dual,primal) obj_ref - f  }; 
z0  = [];   % we don't have a good guess for the dual
tic;
[ x, out, optsOut ] = solver_sBP( A, b, mu, x0, z0, opts );
time_TFOCS = toc;

fprintf('for (original signal, IPM solution, TFOCS solution),\n   NNZ:\n\t%d\t\t%d\t\t%d\n',...
    nnz(x_original),nnz(x_ref), nnz(x) );
fprintf('   error vs. original, rel. l2 norm:\n\t%.2e\t%.2e\t%.2e\n',...
    0, er_signal(x_ref), er_signal(x) );
er_signal1 = @(x) norm(x-x_original,Inf);
fprintf('   error vs. original, lInf norm:\n\t%.2e\t%.2e\t%.2e\n',...
    0, er_signal1(x_ref), er_signal1(x) );
fprintf('   time to solve:\n\tN/A\t\t%.1fs\t\t%.1fs\n',...
    time_IPM, time_TFOCS );

%% Here are some alternative ways to call it
opts = [];
opts.maxIts     = 500;

A_TFOCS = linop_matrix( A );    % not necessary, but one way to do it
[ x, out, optsOut ] = solver_sBP( A_TFOCS, b, mu, x0, z0, opts );

% We can also pass in function handles
Af  = @(x) A*x;
At  = @(y) A'*y;
A_TFOCS = linop_handles( [M,N], Af, At );
[ x, out, optsOut ] = solver_sBP( A_TFOCS, b, mu, x0, z0, opts );

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
