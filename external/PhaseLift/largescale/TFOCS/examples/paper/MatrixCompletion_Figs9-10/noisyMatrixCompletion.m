%{

Like noiselessMatrixCompletion except it includes noise

Recreates fig 10 in the paper

Unlike the noiseless version, this is very fast to run, since we have
very small matrices (that was a necessity, since we needed to run an IPM
to compute the reference solution).

It should take about 10 seconds to run this whole file

%}
%% for speedier computation, use PROPACK

if ~license('test','communication_toolbox')
    awgn = @(x,snr,ignore) x + ...
        10^( (10*log10(sum(abs(x(:)).^2)/length(x(:))) - snr)/20 )*randn(size(x));
end
%%

randn('state',sum('NuclearNorm'));
rand('state',sum('NuclearNorm2'));

% hard = false;
% SNR = Inf;
% oversample = 5;

% The setting used to make the figure in the paper:
M = 50; N = 45; R = 20; hard=true; SNR = 30; oversample = 1.3;% aka "hardcase"

df = R*(M+N-R);
Left  = randn(M,R);
Right = randn(N,R);
k = round(oversample*df); k = min( k, round(.8*M*N) );
omega = randperm(M*N);
omega = sort(omega(1:k)).';
p = k/(M*N);
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

b = b_exact;
if SNR < Inf
    disp('Adding noise...');
    b = awgn(b,SNR,'measured');
    EPS = norm(b-b_exact);
else
    EPS = 0;
end
nuclear_norm = @(X) sum(svd(X,'econ'));

if M*N > 50*50
    R = input('Matrix is large -- you sure you want to run CVX?','s');
    if strcmpi(s,'no') || strcmpi(s,'n')
        error('matrix too large');
    end
end

if M==50 && N==45 && R==20 && hard==true && SNR == 30
    disp('Loading pre-computed solution from file');
    load cvxSolution_hard_Aug26
    CVX = true;
    normXcvx = norm(Xcvx,'fro');
    err = @(X) norm( X - Xcvx, 'fro') / normXcvx;
    fStarBase   = nuclear_norm( Xcvx );
elseif true % (change this to "false" if you don't want a reference soln)
    % compute a solution via CVX
    
    % Note: CVX, on a 80 x 90, rank 50, problem, took 2.95 hours
    % For 50 x 45 matrix, 1800 entries, took 5.4 minutes
    tic
    CVX = true;
    if EPS
      cvx_begin
        variable Xcvx(M,N)
        minimize norm_nuc(Xcvx)
        subject to
            norm(Xcvx(omega) - b ) <= EPS
	  cvx_end
    else
      cvx_begin
        variable Xcvx(M,N)
        minimize norm_nuc(Xcvx)
        subject to
            Xcvx(omega) == b
	  cvx_end
    end
    fprintf('Error with cvx is %.2e\n', er1(Xcvx) );
    toc
    save cvxSolution_hard_Aug26 Xcvx omega R hard b_exact b Left Right
    normXcvx = norm(Xcvx,'fro');
    err = @(X) norm( X - Xcvx, 'fro') / normXcvx;
    fStarBase   = nuclear_norm( Xcvx );
    
else
    CVX = false;
    Xcvx = [];
    normXcvx = [];
    
    % we assume the low-rank matrix is the optimal matrix...
    % (for noisy measurements, this is very unlikely).
    err     = @(X) norm( X - X_exact,'fro')/normX_exact;
    fStarBase   = nuclear_norm(X_exact );
end

%% set some parameters

L = 1;  % bound on Lipschitz constant

optsBase = [];
optsBase.maxits     = 500;
if hard,
    optsBase.maxits = 300;
end
optsBase.maxmin     = -1;  % tell it to do maximization of dual
optsBase.saddle     = 1;
optsBase.printEvery = 20;
optsBase.L0         = L;
optsBase.tol        = 0;
optsBase.continuation = false;
optsBase.stopCrit   = 2; % legacy stopping criteria


%--  Setup tests --
testOptions = {};

% opts = optsBase;
% opts.name   = 'GRA, t=1/L';
% opts.L0     = L;
% opts.Lexact = L;
% opts.beta = 1;
% opts.alg = 'GRA';
% testOptions{end+1} = opts;
% 
% opts = optsBase;
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

% opts = optsBase;
% opts.name   = 'AT, t=1/L';
% opts.L0     = L;
% opts.Lexact = L;
% opts.beta = 1;
% opts.alg = 'AT';
% testOptions{end+1} = opts;

opts = optsBase;
opts.name   = 'AT, backtracking';
opts.L0      = L;
opts.alg = 'AT';
testOptions{end+1} = opts;

% opts = optsBase;
% opts.name   = 'AT, restart every 10';
% opts.L0      = L;
% opts.alg    = 'AT';
% opts.restart = 10;
% testOptions{end+1} = opts;

% opts = optsBase;
% opts.name   = 'AT, restart every 50';
% opts.L0      = L;
% opts.alg = 'AT';
% opts.restart = 50;
% testOptions{end+1} = opts;
% 
% opts = optsBase;
% opts.name   = 'AT, restart every 100';
% opts.L0      = L;
% opts.alg = 'AT';
% opts.restart = 100;
% testOptions{end+1} = opts;

opts = optsBase;
opts.name   = 'AT, continuation';
opts.L0      = L;
opts.alg = 'AT';
opts.continuation = true;
testOptions{end+1} = opts;

%% Solve

testErrors      = {};
times           = [];
shrinkageCalls  = [];
continuationLocations = {};
Y = sparse(omegaI,omegaJ,b,M,N); % be careful: this is modified in place
% kicking:
nY = normest(Y);



for test_i = 1:length( testOptions )

    mu = .0005;

    opts = testOptions{test_i};
    
    CONTINUATION = opts.continuation;
    opts = rmfield(opts,'continuation');
    ERR = [];
    k_end = 1;
    if isfield(opts,'restart') 
        k_end = round( opts.maxits / opts.restart );
        opts.maxits = opts.restart;
    elseif CONTINUATION
        k_end   = 10;
        mu      = .01;
        opts.tol= 1;
    end
   
    opts.restart   = Inf;
    lambda0 = b/nY; 
    lambda = lambda0;
    Xplug   = 0;
%     Xplug   = Y;  % this improves performance a bit...
    Xold    = 0;
    
    if isfield(opts,'L0') && ~isempty(opts.L0)
        opts.L0     = opts.L0/mu;  % legacy code was implicitly scaled
    end
    if isfield(opts,'Lexact') && ~isempty(opts.Lexact)
        opts.Lexact     = opts.Lexact/mu;  % legacy code was implicitly scaled
    end

    prox_nuclear('reset'); sc = 0;
    contLoc = [];

    tic
    for k = 1:k_end
        
        if CONTINUATION
            opts.tol = opts.tol/2;
            if k == k_end
                opts.tol = min(1e-3, opts.tol );
            end
        end
        
        opts.errFcn     = {@(f,dual,x)err(x)};
        [X, out, optsOut ] = solver_sNuclearBPDN({M,N,omega}, b,EPS, mu,Xplug,lambda, rmfield(opts,'name') );
        lambda = out.dual;
        sc_k = prox_nuclear('reset');
        sc = sc + sc_k;
        
        if CONTINUATION, 
            Xplug = X; 
%             Xplug = X + (k-1)/(k+2)*(X-Xold);  Xold = X;  % accelerated
        end
        
        ERR = [ERR; out.err ];
        contLoc = [contLoc, out.niter];
        
    end
    t=toc
    testErrors{end+1} = ERR;
    times   = [times;t];
    shrinkageCalls = [shrinkageCalls; sc ]; % a more fair comparison than just the # of iterations
    continuationLocations{end+1} = contLoc;

end
%% sample plot
figure(1);
semilogy( ERR(:,1) )

hold all

%% nice plot with all the experiments
% We can choose to plot vs. # iterations or vs. # shrinkage calls (aka
% SoftThresholdSingVal )
plotIter = false;  % whether to plot vs. # iterations or not

figure(); co = get(gca,'colororder'); close;
co1 = co([3:5],:);

figure(2);
clf;
set(gcf,'DefaultAxesColorOrder', circshift(co1,-1) );
handles = []; rates = [];
legendString = {};
iter = 1:optsBase.maxits;
for test_i = 1:length( testOptions )
    er=testErrors{test_i}(:,1);
    nn = length(er);
    cLocs = cumsum(continuationLocations{test_i});
    if plotIter
        xGrid = 1:nn;
    else
        xGrid = linspace(1,shrinkageCalls(test_i),nn);
    end
    h=semilogy( xGrid, er ,'linewidth',2);
    hold all
    legendString{end+1} = testOptions{test_i}.name;
    handles= [handles,h];
    if length(cLocs) > 1 
        clr = get(h,'color');
        semilogy( xGrid(cLocs), er(cLocs),'o','color',clr,'markersize',10 );
    end
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
xlim([0,700]);
ylim( [1e-6,1] );

%%
title(sprintf('M=%d, N=%d, p=%.1f%%, d.o.f. = %.f%%, R=%d, \\mu=%.2e, L=1',M,N,100*p,100*df/(M*N),R,mu) );
