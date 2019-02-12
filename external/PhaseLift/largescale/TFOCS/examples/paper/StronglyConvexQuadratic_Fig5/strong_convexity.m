%% Tests the solvers on a simple unconstrained quadratic function

%{
    This script generates and solves a model of the form
        minimize c' * x + x' * D * x / 2
    where D is positive definite. It first converts it to the form
        minimize ( x - x_star )' * D * ( x - x_star ) / 2 + f_star
    where x_star = -D^{-1}c and f_star = c' * x_star / 2, so that we
    can measure optimality errors very precisely.

    This model is strongly convex with parameter min(eig(D)). For
    that reason, optimal first-order methods perform poorly unless
    tuned to exploit that strong convexity. As the experiment shows,
    restart methods can be tuned to perform approximately as well
    as an optimally tuned N83 method.

    The test should take about 40 seconds

    The new figure isn't identical to the one in the paper because
    we've improved the backtracking a little bit, but there 
    are no qualitative differences.
%}

randn( 'state', sum('quadratic test') );
N = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCT THE TEST PROBLEM %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = randn(N,1);
DIAGONAL = false;
if DIAGONAL,
    D = abs(randn(N,1));
    s = D;
    x_star = - D .\ c;
else
    D = randn(N,N);
    D = D * D';
    s = svd(D);
    x_star = - D \ c;
end
L = max(s); 
mu = min(s);
f_star = 0.5 * c' * x_star;
n_star = norm( x_star );
x0 = zeros(N,1);
smoothF = smooth_quad( D, c );
x_error = @(x)norm(x-x_star)/n_star;
f_error = @(f)max((f-f_star)/abs(f_star),realmin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET UP THE TEST PARAMETERS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Since publication we've added some enhancements to the strong
% convexity support, and these tests reflect that. If you want to
% reproduce the figure from the paper, set recreateFigs = true
recreateFigs = true;

opts = [];
opts.tol        = 0;
opts.maxits     = 3000;
opts.maxmin     = 1;
opts.printEvery = 500;
opts.beta       = 0.5;
opts.restart    = Inf;
opts.errFcn     = { @(f,x) x_error(x), @(f,x) f_error(f) };
opts.L0         = L;
opts.Lexact     = L;

constantStepSize = false;
if constantStepSize
    opts.beta   = 1; % turns off backtracking
    
    % for gradient descent, we can actually use the constant-step-size
    % of 2/(L+mu) (> 1/L) and still have convergence (and of course,
    % this gives us a better rate)
end


% Setup tests
testOptions = {};
legendStrings = {};

opts.alg = 'GRA';
testOptions{end+1} = opts;
legendStrings{end+1} = 'GRA';

opts.alg = 'AT';
testOptions{end+1} = opts;
legendStrings{end+1} = 'AT, no restart';

opts.restart = 5;
testOptions{end+1} = opts;
legendStrings{end+1} = 'AT, restart 5';

opts.restart = 10;
testOptions{end+1} = opts;
legendStrings{end+1} = 'AT, restart 10';

opts.restart = 50;
testOptions{end+1} = opts;
legendStrings{end+1} = 'AT, restart 50';

opts.restart = 100;
testOptions{end+1} = opts;
legendStrings{end+1} = 'AT, restart 100';

if ~recreateFigs
    
    % Since publication we've added some enhancements to the strong
    % convexity support, and these tests reflect that.
    
    opts.restart = -Inf;
    testOptions{end+1} = opts;
    legendStrings{end+1} = 'AT, no regress';

    opts.L0      = 1;
    opts.Lexact  = L * 1e6;
    opts.mu      = mu * 1e6;
    opts.restart = Inf;
    testOptions{end+1} = opts;
    legendStrings{end+1} = 'AT, using m';

else
    
    % Select this branch to see the original results.
    
    opts.beta    = 1;  % turn off backtracking
    opts.Lexact  = L;
    opts.mu      = mu;
    opts.restart = Inf;
    testOptions{end+1} = opts;
    legendStrings{end+1} = 'N83, using m';

end

%%%%%%%%%%%%%%%%%
% RUN THE TESTS %
%%%%%%%%%%%%%%%%%
tic
errors = zeros(opts.maxits,2,length(testOptions));
for test_i = 1 : length(testOptions),
    opts = testOptions{test_i};
    fprintf( 'Test: %s\n', legendStrings{test_i} );
    [ x, out, optsOut ] = tfocs( smoothF, [], [], zeros(N,1), opts );
    errors(1:size(out.err,1),:,test_i) = out.err;
end
toc

%%%%%%%%%%%%%%%%%%%%
% PLOT THE RESULTS %
%%%%%%%%%%%%%%%%%%%%
figure();
co = get(gca,'colororder');
close;

figure(1); clf;
set(gcf,'DefaultAxesColorOrder', circshift(co,4) ); % make colors agree with the plot in the paper

h=semilogy( squeeze(errors(:,1,:)),'linewidth',2 );
xlabel('iterations $k$','interpreter','latex','fontsize',16);
ylabel('$\|x_k\textrm{--}x^\star\|/\|x^\star\|$','interpreter','latex','fontsize',16);
hl=legend( legendStrings{:}, 'location','southeast');
set(hl,'fontsize',20);
ylim( [1e-13,1] );
title('Error $||x_k \mathrm{--} x^\star||$','interpreter','latex');
gra_indx = find( strcmpi(legendStrings,'GRA') );
set(h(gra_indx),'linestyle','--');

figure(2);
set(gcf,'DefaultAxesColorOrder', circshift(co,4) ); % make colors agree with the plot in the paper
h=semilogy( squeeze(errors(:,2,:)),'linewidth',2 );
xlabel('$k$','interpreter','latex');
ylabel('$(f_k\mathrm{--}f^\star)/f^\star$','interpreter','latex');
legend( legendStrings{:}, 'location','northeast');
ylim( [1e-13,1] );
title('Error $f_k \mathrm{--} f^\star$','interpreter','latex');
gra_indx = find( strcmpi(legendStrings,'GRA') );
set(h(gra_indx),'linestyle','--');

%% print to PDF
figure(1);
a = -.2;
orient landscape
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf,'PaperPosition',[a a (11-a) (8.5-a)]);
print('-dpdf','fig5_recreation');
