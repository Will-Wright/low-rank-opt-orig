%{
    Demonstrates the exact penalty property of the Dantzig selector.
    Constructs an unsmoothed DS model with a known exact solution, then
    solves a sequence of smoothed problems with a range of values of mu.
    This generates figure 2 from the paper.

    There are a lot of points, and each one is solved to high precision,
    so this test may take a while -- i.e. 8 minutes or so.
%}

if exist('dtest_10dB.mat');
    load dtest_10dB.mat
    SS.type = '()'; SS.subs = { omega, ':' };
    N = length(xExact); M = length(b);
    Af = @(x) subsref(dct(x),SS);
    At = @(x) idct(subsasgn(zeros(N,size(x,2)),SS,x));
    xPlug = xPlug_exact;
else
    opts          = [];
    opts.seed     = 23421;
    opts.smoothed = false;
    opts.SNR      = 10; 
    opts.type     = 'dct';
    opts.N        = 2^11;
    opts.M        = round(opts.N/4);
    [xExact,A,b,normA2,delta0,D,x0,xPlug,mu0,omega] = setup_Dantzig(opts);
    N = length(xExact); M = length(b);
end
A = linop_handles( [M,N], Af, At );

opts            = [];
opts.tol        = 1e-11;
opts.maxits     = 5000;
opts.printEvery = 1000;
opts.restart    = 400;
opts.errFcn     = @(f,dual,primal) norm(primal-xExact)/norm(xExact);
z0              = [];
x0              = [];

muList = logspace( -1.49, -1.475, 4 );
muList = [ muList, logspace( -1.52, -1.3, 10 ) ];
muList = sort( [ muList, logspace( -2, 0, 10 ) ], 'descend' );
errList = zeros(size(muList));
fprintf('Beginning test at %s\n', datestr(now) );
tic
for mu_i = 1:length(muList)
    mu = muList( mu_i );
    fprintf( 'Solving with mu = %.2e\n', mu );
    [ x2, out ] = solver_sDantzig( { A, D }, b, delta0, mu, x0, z0, opts );
    errList(mu_i) = norm(x2-xExact)/norm(xExact);
    fprintf( 'Error: %.2e\n', errList(mu_i) );
    z0 = out.dual; % warm-start the next iteration (but do NOT update x0)
end
fprintf('Finished test at %s\n', datestr(now) );
toc

clf
loglog( muList, errList, 'o-', 'markersize', 10, 'markerfacecolor', 'b' );
xlabel( '$\mu$', 'fontsize', 16, 'interpreter', 'latex' );
ylabel( '$\|x_\mu-x^\star\|/\|x^\star\|$', 'fontsize', 16, 'interpreter', 'latex' );
% print -depsc2 exact_penalty.eps
