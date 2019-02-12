%{
    Tests nuclear norm minimization,

    min_X ||X||_*
s.t.
    ||P(X) - b||_2 <= eps  (or P(X) == b)

where X is a M x N matrix, and P is an observation operator

The solvers solve a regularized version, using
||X||_* + mu/2*||X-X_0||_F^2
and you can either decrease the value of mu,
or use continuation (e.g. accelerated continuation)
to get the solution to the un-regularized problem.

Note: accelerated continuation works very well for matrix completion because with mu
large, the iterates stay very low-rank.
(i.e. updating X0, keeping mu fixed, is a good idea.  Decreasing mu
  is a bad idea).


Reference solution:
    Because this is a noiseless test, if we are sufficiently low-rank
    and have enough observations, then the original low-rank matrix (which
    is known) is also the minimum nuclear norm solution.
    We can't guarantee that we are sufficiently low-rank, but the empirical
    evidence that this is true is overwhelming.

Recreates figure 9 in the paper

Takes about 30 minutes to run this experiment

If you add the path to PROPACK, then it will use PROPACK for the SVD computation,
otherwise it will use MATLAB's builtin "SVDS".
PROPACK is free software by Rasmus Munk Larsen
http://soi.stanford.edu/~rmunk/PROPACK/

Stephen Becker, July 2010, srbecker@caltech.edu

%}

%% for speedier computation, use PROPACK
addpath PROPACK/

%%

randn('state',sum('NuclearNorm'));
rand('state',sum('NuclearNorm2'));

M = 1000; N = 1000; R = 10;  % used in the paper (Fig. 9)
% M = 100; N = 110; R = 2;  % for testing
% M = 300; N = 310; R = 5;  % for testing

df      = R*(M+N-R);
oversample = 5;
Left    = randn(M,R);
Right   = randn(N,R);
k       = round(oversample*df); k = min( k, round(.8*M*N) );
omega   = randperm(M*N);
omega   = sort(omega(1:k)).';
p       = k/(M*N);
fprintf('%d x %d rank %d matrix, observe %d = %.1f x df = %.1f%% entries\n',...
    M,N,R,k,k/df,p*100);
[omegaI,omegaJ] = ind2sub([M,N],omega);
mat = @(x) reshape(x,M,N);
vec = @(x) x(:);
if ~issorted(omega)
    error('Algorithm assumes omega is sorted, for simplicity');
end

% If we were testing a huge matrix, we would only want to check
%   our error on a subset since otherwise the error computation
%   would be too slow.  But for this size matrix, it's not a big deal.
omega_test = 1:(M*N);
omega_test = omega_test(:); % make it a column vector


X_exact = Left*Right';
b_exact = X_exact(omega);  % no noise for now
X_exact_test = X_exact( omega_test );
normX_exact_test = norm(X_exact_test);
normX_exact = norm(X_exact,'fro');

EPS = 0;
b = b_exact;  % no noise

nuclear_norm = @(X) sum(svd(X,'econ'));

if M*N < 50*50
    CVX = true;
    % get reference via CVX
      cvx_begin
        variable Xcvx(M,N)
        minimize norm_nuc(Xcvx) %+ mu/2*sum_square( Xcvx(:) - Xplug(:) )
        subject to
            Xcvx(omega) == b
	  cvx_end
    fprintf('Error with cvx is %.2e\n', norm(Xcvx-X_exact,'fro')/norm(X_exact,'fro'));
    normXcvx = norm(Xcvx,'fro');
    err = @(X) norm( X - Xcvx, 'fro') / normXcvx;
    fStarBase   = nuclear_norm( Xcvx );
else
    CVX     = false;
    Xcvx    = [];
    normXcvx = [];
    err     = @(X) norm( X - X_exact,'fro')/normX_exact;
    fStarBase   = nuclear_norm(X_exact );
end

%% set some parameters

Xplug   = zeros(M,N); normXplug = 0;
L       = 1;  % bound on Lipschitz constant

optsBase = [];
optsBase.maxits     = 500;      % Fig 9 in the paper used 500
% optsBase.maxits     = 25; disp('Warning: using very low maxits'); % for testing
optsBase.maxmin     = -1;        % tell it to do maximization of dual
optsBase.saddle     = 1;
optsBase.printEvery = 10;
% optsBase.dualAvg    = false;  % since we use packed storage
optsBase.L0         = L;
% optsBase.Lexact     = L;
% optsBase.L          = p/1.2; % SVT convention
optsBase.alg     = 'AT';
optsBase.tol        = 0;


%--  Setup tests --
% To add a test to the list, just uncomment out the relevant section

testOptions = {};

% opts = optsBase; % this one is not used in Fig 9 in the paper
% opts.name   = 'GRA, t=1/L';
% opts.Lexact = L;
% opts.beta = 1;
% opts.alg = 'GRA';
% testOptions{end+1} = opts;

% opts = optsBase; % this one is not used in Fig 9 in the paper
% opts.name   = 'GRA, t=p/1.2';
% opts.L0     = p/1.2;
% opts.Lexact = p/1.2;
% opts.beta = 1;
% opts.alg = 'GRA';
% testOptions{end+1} = opts;


opts = optsBase;
opts.name   = 'GRA, backtracking';
opts.L0     = L;
opts.alg = 'GRA';
testOptions{end+1} = opts;

% opts = optsBase; % this one is not used in Fig 9 in the paper
% opts.name   = 'AT, t=1/L';
% opts.Lexact = L;
% opts.beta = 1;
% opts.alg = 'AT';
% testOptions{end+1} = opts;

opts = optsBase;
opts.name   = 'AT, backtracking';
opts.alg = 'AT';
testOptions{end+1} = opts;

% opts = optsBase;      % this one is not used in Fig 9 in the paper
% opts.name   = 'AT, restart every 5';
% opts.alg = 'AT';
% opts.restart = 5;
% testOptions{end+1} = opts;

opts = optsBase;
opts.name   = 'AT, restart every 10';
opts.alg = 'AT';
opts.restart = 10;
testOptions{end+1} = opts;

opts = optsBase;
opts.name   = 'AT, restart every 50';
opts.alg = 'AT';
opts.restart = 50;
testOptions{end+1} = opts;

opts = optsBase;
opts.name   = 'AT, restart every 100';
opts.alg = 'AT';
opts.restart = 100;
testOptions{end+1} = opts;

%% Solve
mu = 1e-4; % use mu = 1e-4 in Fig. 9 in the paper

fStar = fStarBase + mu/2*norm(X_exact-Xplug,'fro')^2; 
testErrors      = {};
times           = [];
shrinkageCalls  = [];
lambda0     = [];

for test_i = 1:length( testOptions )
% for test_i = 1:1
    opts = testOptions{test_i};
    
    opts.errFcn     = {@(f,dual,x)err(x), @(f,dual,x) fStar - f};

        
    % Old-style restart (in newer code, it is builtin)
    lambda = lambda0;
    k_end = 1;
    ERR = [];
    if isfield(opts,'restart')
        k_end = round( opts.maxits / opts.restart );
        opts.maxits = opts.restart;
    end
    opts.restart = Inf; % override the defaults
    
    if isfield(opts,'L0') && ~isempty(opts.L0)
        opts.L0     = opts.L0/mu;  % legacy code was implicitly scaled
    end
    if isfield(opts,'Lexact') && ~isempty(opts.Lexact)
        opts.Lexact     = opts.Lexact/mu;  % legacy code was implicitly scaled
    end
    
    
    % for debugging:
    %     opts.printEvery = 1; opts.maxits = 10;
    %     opts.beta       = 1; % turn off backtracking
    %     opts.errFcn{end+1} = @(f,dual,x) rank(x);  % SLOW! Comment this out

    prox_nuclear('reset'); sc = 0;
    tic
    for k = 1:k_end
        [X, out, optsOut ] = solver_sNuclearBP({M,N,omega}, b, mu,Xplug,lambda, rmfield(opts,'name') );
        ERR = [ERR; out.err ];
        lambda = out.dual;
        sc_k = prox_nuclear('reset');
        sc = sc + sc_k;
    end
    t=toc
    testErrors{end+1} = ERR;
%     testErrors{end+1} = out.err;
    times   = [times;t];
    
    shrinkageCalls = [shrinkageCalls; sc ]; % a more fair comparison than just the # of iterations

end
%% sample plot
% figure(1);
% semilogy( out.err(:,1) )
% 
% hold all

%% nice plot with all the experiments
% We can choose to plot vs. # iterations or vs. # shrinkage calls (aka
% SoftThresholdSingVal )
plotIter = false;  % whether to plot vs. # iterations or not


figure(); co = get(gca,'colororder'); close;
co1 = co([1:2,4:6],:);
co1(end,:) = [0,0,0];
figure(2);
clf;
set(gcf,'DefaultAxesColorOrder', circshift(co1,-2) );
handles = []; rates = [];
legendString = {};
iter = 1:optsBase.maxits;
for test_i = 1:length( testOptions )
    er=testErrors{test_i}(:,1);

    nn = length(er);
    if plotIter
        xGrid = 1:nn;
    else
        xGrid = linspace(1,shrinkageCalls(test_i),nn);
    end
    h=semilogy( xGrid, er ,'linewidth',2);
    hold all
    legendString{end+1} = testOptions{test_i}.name;
    handles= [handles,h];
    if ~isempty( strfind(legendString{end},'GRA') )
        set(h,'linestyle','--');
    end
    
end

legend( handles, legendString,'fontsize',20 ,'location','southwest');

if plotIter
    xlabel('iterations','fontsize',16);
else
    xlabel('calls to SoftThresholdSingVal','fontsize',16);
end
ylabel('error','fontsize',16);
xlim([0,1200])
%%
title(sprintf('M=%d, N=%d, p=%.1f%%, d.o.f. = %.f%%, R=%d, \\mu=%.2e, L=1',M,N,100*p,100*df/(M*N),R,mu) );
