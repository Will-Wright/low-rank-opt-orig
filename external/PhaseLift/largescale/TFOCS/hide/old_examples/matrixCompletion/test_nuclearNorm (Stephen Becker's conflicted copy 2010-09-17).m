%{
    Tests nuclear norm minimization,

    min_X ||X||_*
s.t.
    P(X) == b, or || P(X) - b || <= eps

where X is a M x N matrix, and P is an observation operator

The solvers solve a regularized version, using
||X||_* + mu/2*||X-X_0||_F^2


TODO: generate problem with CVX and save as a .mat file, so the user doesn't need CVX

This still needs to be cleaned up!

%}
nuclear_norm = @(x) sum(svd(x,'econ'));

% Try to load the problem from disk
fileName = 'nuclearNorm_problem1_noiseless';
if exist([fileName,'.mat'],'file')
    load(fileName);
    fprintf('Loaded problem from %s\n', fileName );
else
    
    % Generate a new problem

    randn('state',sum('NuclearNorm'));
    rand('state',sum('NuclearNorm2'));
    
    M = 30; N = 30; R = 2;

    df = R*(M+N-R);
    oversample = 5;
    Left  = randn(M,R);
    Right = randn(N,R);
    k = round(oversample*df); 
    k = min( k, round(.8*M*N) );
    omega = randperm(M*N);
    omega = sort(omega(1:k)).';

    X_original = Left*Right';       % the "original" signal -- may not be optimal value
    b_original = X_original(omega); 
    EPS = 0;        % noise level
    b = b_original;
    obj_original = nuclear_norm(X_original);

    % get reference via CVX
    cvx_begin
        cvx_precision best
        variable Xcvx(M,N)
        minimize norm_nuc(Xcvx)
        subject to
            Xcvx(omega) == b
    cvx_end
    X_reference = Xcvx;         % the nuclear norm minimizer
    obj_reference = nuclear_norm(X_reference);
    
    save(fileName,'X_original','X_reference','omega','b','obj_original',...
        'Left','Right','EPS','b_original','R','obj_reference');
    fprintf('Saved data to file %s\n', fileName);
    
end

[M,N]           = size(X_reference);
norm_X_reference = norm(X_reference,'fro');
er_reference    = @(x) norm(x-X_reference,'fro')/norm_X_reference;
resid           = @(x) norm(x(omega)-b)/norm(b);  % change if b is noisy

[omegaI,omegaJ] = ind2sub([M,N],omega);
mat = @(x) reshape(x,M,N);
vec = @(x) x(:);
    
k  = length(omega);
p  = k/(M*N);
df = R*(M+N-R);
fprintf('%d x %d rank %d matrix, observe %d = %.1f x df = %.1f%% entries\n',...
    M,N,R,k,k/df,p*100);
fprintf(' Nuclear norm solution and original matrix differ by %.2e\n',...
    norm(X_reference-X_original,'fro')/norm_X_reference );

%% set some parameters

L = 1;  % bound on Lipschitz constant
optsBase = [];
optsBase.maxits     = 500;      
optsBase.maxmin     = -1;   % tell it to do maximization of dual
optsBase.saddle     = 1;
optsBase.printEvery = 10;
optsBase.L          = L;
optsBase.solver     = @solver_AT;
optsBase.tol        = 0;

testOptions = {};
opts = optsBase;
opts.name   = 'GRA, backtracking';
opts.L      = L;
opts.solver = @solver_GradientDescent;
testOptions{end+1} = opts;

%%
mu      = .001;
Xplug   = zeros(M,N); normXplug = 0;
normXplug       = norm(Xplug,'fro')^2;
fStar           = obj_reference + mu/2*norm(X_reference-Xplug,'fro')^2;

Y = sparse(omegaI,omegaJ,b,M,N); % be careful: this is modified in place
% kicking:
nY = my_normest(Y,1e-4,300);
lambda0 = b/mu/nY;
%% Run it, very simply:
% addpath ~/Dropbox/TFOCS/
[ x, z, out, opts ] = solver_sNuclearBP( {M,N,omega}, b, mu, Xplug, [] ); 





%% For now, ignore everything below here

%%
testErrors      = {};
times           = [];
shrinkageCalls  = [];


for test_i = 1:length( testOptions )
    opts = testOptions{test_i};
    solver = opts.solver;   % Change this
    lambda = lambda0;
    ERR = [];
    k_end = 1;
    if isfield(opts,'restart') 
        k_end = round( opts.maxits / opts.restart );
        opts.maxits = opts.restart;  % Change this
    end
    shrinkMatrix();  % zero it out.   Change this
    tic
    for k = 1:k_end
        
%         opts.errFcn        = { @(f,l,x) norm(unpackSVD(x)-Xcvx)/normXcvx };
        opts.errFcn     = { @(f,l,x) er_reference(x), @(f,l,x) fStar - f };
   
        % Old:
        observations = @(X) X(omega);
        opts.errFcn{end+1} = @(l,f,x) (norm( observations(unpackSVD(x)) - b )- EPS)/norm(b); % residual
        
   
        
        % this operates on a packed matrix
        if norm(normXplug,'fro') < 10*eps
            obj = @(x) sum(unpackSVD(x,'s')) + mu/2*norm(unpackSVD(x,'s'))^2;
            opts.errFcn{end+1} = @(l,f,x) obj(x) - fStar;
            % smooth_nuclear computes it more efficiently
        end
        
        rank = @(x) unpackSVD(x,'rank');
        opts.errFcn{end+1} = @(l,f,x) rank(x);
        
     
        % For both, use this:
        %Y = sparse(omegaI,omegaJ,b,M,N); % be careful: this is modified in place
        %linear = [];
        %smooth = @(varargin) smooth_nuclear(Y,b,mu,XplugF, varargin{:} );
        %proj   = @(varargin) nonsmooth_nuclear(EPS,varargin{:}); % Psi.  Same as for BPDN!
        %
        %
        %[lambda,out,optsOut] = solver( smooth, linear, proj, lambda, opts );
%
        %x = out.dual;
        %[U,S,V] = unpackSVD(x);
        %s = diag(S);
        %X = U*S*V';

            % Q: does optsOut still work?
        [x,lambda,out] = solver_sNuclearBP( omega, b, my, Xplug, [], opts );



        %     Xplug = X;
        %     XplugF = Xplug;
        % Make Xplug with an implicit form:
        %     XplugF = { @(x) U*(S*(V'*x)), @(x) V*(S*(U'*x)), norm(s)^2 }; normXplug = norm(s)^2;
        %     XplugF = X + (k-1)/(k+2)*( X - Xplug );normXplug=norm(XplugF,'fro');
        % will this work if we keep it as low-rank multiply?
%         lambda0 = lambda;
        %     opts.tol = opts.tol/2;
        
        ERR = [ERR; out.err ];
        %     dmu = 1;
        %     mu = dmu*mu;
        %     lambda0 = lambda0/(dmu);
        
    end
    t=toc
    testErrors{end+1} = ERR;
    times   = [times;t];
    sc      = shrinkMatrix();
    shrinkageCalls = [shrinkageCalls; sc ]; % a more fair comparison than just the # of iterations

end
%% sample plot
figure(1);
semilogy( ERR(:,1) )

hold all

%% nice plot with all the experiments
% We can choose to plot vs. # iterations or vs. # shrinkage calls (aka
% SoftThresholdSingVal )
plotIter = false;  % whether to plot vs. # iterations or not

% can also choose whether to plot the estimated error rates
plotRates = false;


figure(2);
clf;
handles = []; rates = [];
legendString = {};
iter = 1:optsBase.maxits;
for test_i = 1:length( testOptions )
    er=testErrors{test_i}(:,1);
    er_rank=testErrors{test_i}(:,3); % for rank
    if any(er_rank>R), disp('found ranks > true rank'); end
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
    
%     r_est = 1-mean( er(2:end)./er(1:end-1) );
    er = er( er > 1e-13) ;
    er = er( er < .5 );
    if ~isempty(er) && er(1) ~= 0
        r_est = 1 - (er(end)/er(1))^(1/length(er));
    else
        r_est = 1;
    end
    rates = [rates;r_est];
    
    if plotRates
        gtext( sprintf('%.4f',r_est), 'fontsize', 16 );
    end
end

legend( handles, legendString,'fontsize',20 ,'location','southwest');

if plotIter
    xlabel('iterations','fontsize',16);
else
    xlabel('calls to SoftThresholdSingVal','fontsize',16);
end
ylabel('error','fontsize',16);
%%
title(sprintf('M=%d, N=%d, p=%.1f%%, d.o.f. = %.f%%, R=%d, \\mu=%.2e, L=1',M,N,100*p,100*df/(M*N),R,mu) );
