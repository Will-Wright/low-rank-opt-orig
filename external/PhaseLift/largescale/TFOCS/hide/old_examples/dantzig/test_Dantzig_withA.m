%{
    Run the first part of "old/test_Dantzig.m" to setup the problem
    and the reference solution.

    This shows how to use Dantzig with the simpler formulation
    (linear function is the identity)
    (this is less efficient)

    Stephen Becker, 2/24/10

Modified:  Wed, Apr 14 2010 Stephen Becker

Funny stuff: L exceeds L_exact sometimes (by a LOT!).  Due to numerical error?

%}

lambda0 = zeros(N,1);
opts.Lexact = norm(diag(D)*A'*A)^2 / mu;
% opts.maxcounts  = [Inf,Inf,10000,Inf];
opts.maxcounts  = [Inf,Inf,Inf,Inf];
opts.maxits     = 5000;
opts.printEvery = 500;      % how often to print to screen
errob = @(f,l,xk) abs( f - oExact ) / oExact;  % requires knowledge of oExact
errd  = @(f,l,xk) norm( l - lambdaExact, inf ) / norm( lambdaExact,Inf );
errp  = @(f,l,xk) norm( xk - xExact, inf ) / norm( xExact,Inf );
errsp = @(f,l,xk) nnz( ~(abs(xk)>1e-8*norm(xk,Inf)) ~= ~xExact );
errsd = @(f,l,xk) nnz( ~(abs(l)>1e-8*norm(l,Inf)) ~= ~lambdaExact );
opts.errFcn = { errob, errd, errp, errsp, errsd };

ax = [+Inf,-Inf,+Inf,-Inf];
% solver = @solver_AT;

%% First, solve the efficient way (as in test_Dantzig.m)
solver = @solver_AT;
% solver = @solver_GradientDescent; opts.L = norm(diag(D)*A'*A)^2 / mu;
opts.printEvery = 500;

odef = opts;
% linear = @(varargin) linear_Dantzig( Af, At, b_orig, D, varargin{:} );
linear = @(varargin) linear_Dantzig( Af, At, b, D, varargin{:} ); % for cleaned-up solution
smo    = @(varargin) smooth_Dantzig( mu, mu * xPlug, varargin{:} );
proj   = @(varargin) nonsmooth_Dantzig( delta0, varargin{:} );
f      = @(varargin) composite( smo, linear, varargin{:} );

[lambda,out,opts] = solver( smo, linear, proj, lambda0, odef );
xk = out.dual;
fprintf('--Solved smoothed problem via first-order method, mu = %.1e --\n',mu);
fprintf( 'Objective (exact-primal,exact,exact-dual): (%.8e,%.8e,%.8e)\n', ...
    oExact-obj(xk), oExact, oExact-out.f(end) );
fprintf( 'l_inf error (primal,dual): (%.2e,%.2e)\n', out.err(end,2), out.err(end,3) );
fprintf( 'support error (primal,dual): (%d,%d)\n', out.err(end,4), out.err(end,5) );
fprintf('\t# of calls to A and At: %d. # calls per iteration: %.2f\n',...
    out.counts(end,3), out.counts(end,3)/out.niter );
outEfficient = out;
%% Now, the simpler way
solver = @solver_AT; opts.L = [];
% solver = @solver_GradientDescent; opts.L = norm(diag(D)*A'*A)^2 / mu;
odef = opts;

linear = [];
smo    = @(varargin) smooth_Dantzig_withA( A,A',D.*(A'*b),D,mu,xPlug, varargin{:} );

[lambda,out,opts] = solver( smo, linear, proj, lambda0, odef );
xk = out.dual;
fprintf('--Solved smoothed problem via first-order method, mu = %.1e --\n',mu);
fprintf( 'Objective (exact-primal,exact,exact-dual): (%.8e,%.8e,%.8e)\n', ...
    oExact-obj(xk), oExact, oExact-out.f(end) );
fprintf( 'l_inf error (primal,dual): (%.2e,%.2e)\n', out.err(end,2), out.err(end,3) );
fprintf( 'support error (primal,dual): (%d,%d)\n', out.err(end,4), out.err(end,5) );
fprintf('\t# of calls to A and At: %d. # calls per iteration: %.2f\n',...
    out.counts(end,3), out.counts(end,3)/out.niter );
outSimple = out;
%% Plot results
figure(1); clf;
subplot(1,2,1);
errType = 1;
semilogy( outSimple.err(:,errType) );
hold all
semilogy( outEfficient.err(:,errType) );
legend('simple','efficient');

% Q: why the discrepancy? Efficient method stagnates
% title('Why does the "simple" implementation not stagnate, when the "efficient" implementation does?');
% ANSWER: I used b for one, and b_orig for the other. That explains it.
%  So the apparent problem is solved.  You can ignore this section.
% (cleanup before releasing code)
xlabel('iterations');
ylabel('error in objective function');

subplot(1,2,2);
errType = 3;
semilogy( outSimple.err(:,errType) );
hold all
semilogy( outEfficient.err(:,errType) );
legend('simple','efficient');
xlabel('iterations');
ylabel('l_\infty error of dual');

%% estimate error rate
figure(2); clf;
errType = 1;
semilogy( outEfficient.err(:,errType) );
hold all
% semilogy( abs(outEfficient.f - oExact)/oExact );
K = 1:outEfficient.niter;

% assuming decay of c/k^2, estimate c by linear regression
e=outEfficient.err(:,errType);
e = e(100:3000);
K = K(100:3000);
c = (K.^-2)*e /sum(K.^-4);
semilogy( K, c./(K.^2) );

% assuming decay of linear rate, estimate rate
e=outEfficient.err(:,errType);K = 1:outEfficient.niter;
e = e(3500:4000);
K = K(3500:4000);
% e = e(100:3000);
% K = K(100:3000);
r = mean( e(2:end)./e(1:end-1) );
a = (r.^K)*e / ( norm(r.^K)^2 );

semilogy( K, a*r.^K )
