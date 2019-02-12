%{
Original: Wed, March 17 2010
Cleaned version: Mon, Oct 25 2010
Stephen Becker, srbecker@caltech.edu

Tests continuation (regular and accelerated),
shows only outer iterations.  So, this isn't the whole story, 
but provides insight.

Test should take about 4 minutes

%}

%%
Create_New_Problem = false;

if Create_New_Problem
    % -- Setup a problem
    randn('state',2010);
    N       = 200;
    M       = 80;
    A       = randn(M,N);
    xExact  = randn(N,1);
    b       = A*xExact;
    epsilon = 1;

    % -- Find a reference solution via CVX (i.e. unregularized)
    cvx_begin
      cvx_quiet(true)
      cvx_precision best
      variable x(N,1)
      minimize norm(x,1)
      subject to
        norm(A*x - b) <= epsilon
    cvx_end
    if ~strcmpi(cvx_status,'Solved')
        fprintf('\n\nError with CVX answer\n\n');
    end
    save ContinuationOuterLoopOnly_experimentForPaper A xExact b epsilon x
else
    % load a problem from disk
    load ContinuationOuterLoopOnly_experimentForPaper
    [M,N] = size(A);
end

x_ref = x;
er =  @(x) norm(x-x_ref)/norm(x_ref);
er2 = @(x) norm(x,1) - norm(x_ref,1);

% -- Repeatedly solve a smoothed problem
mu =  .5; kMax = 30;  % this is what we used in the figure in the paper

t       = clock;
% Either form of c0 works, makes little difference
% c0      = zeros(N,1);
c0      = randn(N,1);  % used this for the paper
for NESTEROV = 0:1
    if NESTEROV
        disp('--- Continuation (accelerated) ---');
    else
        disp('--- Continuation (regular) ---');
    end
    
    c       = c0; 
    x       = zeros(N,1);
    xOld    = x;
    data    = [];
    lambda  = [];
    for k = 1:kMax
        
        if NESTEROV==1 
            c = x + k/(k+3)*(x-xOld);
            xOld = x;
        else
            c = x;
        end
        
        opts        = [];
        opts.tol    = 1e-10;
        opts.maxits = 5000;
        opts.stopCrit = 2; % legacy option in the new TFOCS code
        [x,out]  = solver_sBPDN(A,b,epsilon,mu,c,lambda,opts);
        lambda = out.dual;

        fprintf('Iteration %2d, error is %.2e, ||resid|| is %.2e, ||x||_1-||x^*||_1 is %.2e\n',...
            k, er(x), norm(A*x-b), er2(x) );
        
        % record:
        data(k).er = er(x);
        data(k).er2 = er2(x);
        data(k).resid = norm(A*x-b);
        
    end
    DATA{NESTEROV+1} = data;
end
fprintf('Finished at %s, total wall time is %.2f minutes\n',datestr(now),...
    etime(clock,t)/60 );
%% Plot results
K = 1:kMax;
figure(1); clf;
plotf = @semilogy;

% we can show either data.er ( ||x-x^*|| ) or data.er2 ( f(x) - f(x^*) )

data = DATA{1};
plotf( K, ( [data.er]/data(1).er ).^2,'o-','linewidth',2 )
hold all
data = DATA{2};
plotf( K, ( [data.er]/data(1).er ).^2,'*-','linewidth',2 )

h = legend('Regular continuation (fixed \mu)','Accelerated continuation (fixed \mu)',...\
    'location','southwest');
set(h,'fontsize',22);

if strcmpi( get(gca,'XScale'), 'linear' )
    xlim( [1,kMax] );
else
    xlim( [1,2*kMax] );  % for x on log scale
end

xlabel('outer iteration','fontsize',18)
ylabel('error','fontsize',18)

% print it to a PDF
% orient landscape
% print('-dpdf','AcceleratedContinuation_mu0-5');