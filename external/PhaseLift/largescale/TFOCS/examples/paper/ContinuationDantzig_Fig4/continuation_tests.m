%{

Creates an exact solution to the unsmoothed Dantzig Selector problem,
then solves the smoothed Dantzig Selector problem with various types of continuation.

Figure 4

takes about 7 or 8 minutes to run

%}


% Make sure you ran "setup_path_test" in the parent directory

%% Setup a test problem
opts = [];
opts.seed = 23421;
opts.smoothed = false; 
opts.SNR = 30;
opts.type = 'dct';
opts.N = 2^8;
opts.M = round(opts.N/4);

% fileName = 'continuation_test_OldFormat';
fileName = 'continuation_test_NewFormat';
if ~isempty(strfind(fileName,'New'))
    NewFormat = true;   % using this for the October update of the paper
                        % reference solution is a bit better
else NewFormat = false;
end

if exist([fileName,'.mat'],'file')
    % load file from disk
    fprintf('Loading test problem from %s\n', fileName );
    load(fileName);
else
    fprintf('Cannot find %s.mat, so regenerating test problem\n',fileName);
    % generated the problem again
    if NewFormat
        opts.NewFormat = true;
    end
    [xExact,A,b,normA2,delta0,D,x0,xPlug_exact,mu0,omega,lambda] = setup_Dantzig(opts);
    % This gives us an exact known solution
    save(fileName,'xExact','b','normA2','delta0','D','x0','xPlug_exact',...
        'mu0','omega','lambda' );
end
N = length(xExact); M = length(b);
downsample = @(x) x(omega,:);
SS.type = '()'; SS.subs{1} = omega; SS.subs{2} = ':';
upsample = @(x) subsasgn( zeros(N,size(x,2)),SS,x);
A = {};
A{1} = @(x) downsample(dct(x));
A{2} = @(x) idct(upsample(x));



Aff = A{1}; Att = A{2};

Af = @(x) multiplyA( A{1}, A{2}, 'forward', x );    % this keeps track of # of calls
At = @(x) multiplyA( A{1}, A{2}, 'transpose', x );

% obj = @(x) norm(x,1) + mu0/2*sum_square(x-xPlug_exact);
obj = @(x) norm(x,1);
errp  = @(f,l,xk) norm( xk - xExact) / norm( xExact );

oExact = obj(xExact);
% This error function helps us save error information
globalErrors = @(f,l,xk) ComputeErrors( xExact, oExact, xk, f );

fprintf('  Reference solution violates constraints by  %.2e\n',...
    norm( diag(D)*Att( Aff(xExact) - b),Inf ) - delta0 );
fprintf('  Reference solution has %d nnz, i.e. %.1f%% sparse\n',nnz(xExact), nnz(xExact)/N*100 );
fprintf('  Reference solution has %.1f dB dynamic range\n',...
     20*log10(norm(xExact,inf)/min(abs(xExact(~~xExact))) ) );
%% Prove to the skeptic that we really have a solution to the unsmoothed
%% Dantzig:
%{
Let's verify the KKT conditions.  I will write out the Lagrangian
so that we agree on notation:

primal:  min_x ||x||_1   subject to  || D*A'(b-Ax) ||_inf <= delta0
lagrangian(x,lambda) = ||x||_1 - ( < D*A'(b-Ax), lambda > + <delta0,tau> )
  with constraint ||lambda||_1 <= tau
  (We actually think of tau as exactly ||lambda||_1 )

so, KKT conditions:
  primal feasibility: ||DA'(b-Ax) ||_inf <= delta
  dual feasibiltiy:   ||lambda||_1 <= tau
  complementary slackness: we want a zero primal-dual gap,
    so for the value of the Lagrangian to be equal to the primal objective,
    we need ( < D*A'(b-Ax), lambda > + <delta0,tau> ) = 0
  stationarity:   x must be a stationary point of the Lagragian, i.e.
        0 \in nabla_x lagrangian(x,lambda) , where nabla is the subgradient

We can get this exactly for mu > 0.  For mu = 0, we can still get an 
exact primal solution by taking mu small enough and using the results
of the exact perturbation; however, this doesn't mean that our dual solution
lambda will be valid for mu = 0.  Hence we may see that the final KKT
condition isn't satisfied up to machine precision.

%}
disp('TESTING VIOLATION OF KKT CONDITIONS');
resid = diag(D)*Att( Aff(xExact) - b );
primalInfeas = @(r) norm( r ,Inf ) - delta0 ;
fprintf('KKT 1: primal constraint violated by:\t\t\t%.2e\n', primalInfeas(resid) );
fprintf('KKT 2: ||lambda||_1 = tau.  This is by construction.\t0\n');
primalValue = norm(xExact,1);
compSlackness = resid'*lambda + delta0*norm(lambda,1);
fprintf('KKT 3: complementary slackness (aka duality gap):\t%.2e\n', compSlackness )
supp = @(x) x(find(xExact));
fprintf('KKT 4a: gradient_x(Lagrangian) = 0. On supp(t), norm(gradient):\t\t%.2e\n',...
    norm( supp(   sign(xExact) + Att(Aff(D.*lambda ))  )) );
off_supp = @(x) x(find(~xExact));
fprintf('KKT 4b: On complement of supp(t), ensure subgrad(norm(x,1)) in [-1,1]:\t%.2e\n',...
    norm( off_supp( Att(Aff(D.*lambda ))) , Inf ) - 1 );
%% Clear some variables
ERRS = {};
CONTINUATION_LOCATIONS={};
%% Test the smoothed problem, no continuation
opts=[];
e_l1    = @(xk) abs( norm(xk,1) - norm(xExact,1) );
arg2    = @(list) list(2);
err_l1  = @(f,l,xk)  runningAvg(50,norm(xk,1),'std');
opts.errFcn     = { errp, globalErrors };
stopTol         = 1e-3;
opts.maxits     = 4000;
opts.printEvery = 500;
% opts.solver='solver_AT';  % solver is AT by default

tic
% CONTLIST = 1;
CONTLIST = -2:2;
% CONTLIST = -3:2;
cOffset = 3;
if any(ismember(CONTLIST,-3)), cOffset = 4; end
for CONTINUATION = CONTLIST
  opts.xPlug = zeros(N,1);


  if NewFormat
      mu0     = 1e-2;  %opts.printEvery = 50;
  else
      mu0     = 1e-2;
  end
  tolFinal= 1e-8;
  k_end   = 1;      % for the non-continuation versions
  tol     = tolFinal;
  switch CONTINUATION
    
    case -3
        % no continuation, but with small mu, and using restart
        mu = mu0/10;
    case -2
        % no continuation, but with small mu
        mu = mu0/10;
    case -1
        % no continuation
        mu = mu0;               
    case 0
        % no continuation, but with large mu
        mu = 10*mu0;
    case {1,2}
        % 1: continuation: fixed mu, update x0
        % 2: accelerated continuation, version 1
        tol = norm(b)*1e0;
        k_end = 15;
        mu = 50*mu0;
        mu_continuation = mu;
  end

  xOld            = zeros(N,1);
  opts.tol        = tol;
  opts.lambda0    = zeros(N,1);
  lambdaOld       = opts.lambda0;
  continuation_locations = [];
  rel_differences = [];
  relDiff         = 1;
  error_history   = [];
  multiplyA();  % zero it out.

  for k = 1:k_end
    fprintf('----------------------- Solving with mu = %.1e, k = %d ---------------\n',mu,k);

    stopTol = stopTol/2;
    
    opts.restart = Inf;
    if k == k_end
        opts.tol = tolFinal;
        if CONTINUATION <= 0
            opts.maxits = 10000; % was 12000 before
            opts.restart = Inf;
            if CONTINUATION == -3
                opts.restart = 2000;
            end
        else
            opts.maxits = 5000;
            opts.restart = 1000;  % check this
        end
        opts.stopFcn = [];

    elseif CONTINUATION >= 1
        opts.tol = opts.tol/2;
    end
    if isfield(opts,'stopFcn') && ~isempty(opts.stopFcn)
        opts.tol = 0;
    end
    

    runningAvg(); % zero it out
    A = linop_handles([N,N],Af,At);
    opts.stopCrit = 2;  % legacy stopping criteria
    [ x2, out ] = solver_sDantzig( {A,D}, b, delta0, mu,opts.xPlug,opts.lambda0, opts );
    lambda2 = out.dual;
    
    continuation_locations = [continuation_locations; out.niter ];
    
    switch CONTINUATION
        case 1
            % standard continuation
            opts.lambda0 = lambda2;
            opts.xPlug = x2;
        case 2
            % accelerated continuation
            opts.lambda0 = lambda2;
            opts.xPlug = x2 + (k+0)/(k+0+3)*( x2 - xOld );
    end
    
    relDiff = norm( xOld - x2 )/norm(xOld);
    rel_differences = [rel_differences; relDiff];


    xOld = x2;
    lambdaOld = lambda2;

    error_history = [error_history; out.err ];

  end
  ERRS{CONTINUATION+cOffset} = error_history;
  CONTINUATION_LOCATIONS{CONTINUATION+cOffset} = continuation_locations;

end
toc
%% quick plot
% error_history           = ERRS{1};
% continuation_locations  = CONTINUATION_LOCATIONS{1};
errs = error_history(:,1);
semilogy( 1:length(errs), errs );
%% plot
figure(1);
clf;  
handles = [];
cList   = (-2:2);
for CONTINUATION = cList
    error_history = ERRS{CONTINUATION+cOffset};
    continuation_locations = CONTINUATION_LOCATIONS{CONTINUATION+cOffset};
    errs = error_history(:,1);
    iter = 1:length(errs);
    
    h=semilogy( iter, errs,'linewidth',2 );
    if CONTINUATION <= 0
        set(h,'linestyle','--');
    end
    handles( CONTINUATION + cOffset ) = h;
    hold all
end
hold on

for CONTINUATION = cList
    error_history = ERRS{CONTINUATION+cOffset};
    continuation_locations = CONTINUATION_LOCATIONS{CONTINUATION+cOffset};
    continuation_locations = cumsum(continuation_locations);

    errs = error_history(:,1);
    iter = 1:length(errs);
    if CONTINUATION > 0
        clr = get(handles(CONTINUATION+cOffset),'color');
        semilogy(continuation_locations,errs(continuation_locations),...
            'o','color',clr,'markersize',10 );
    end
end

names = {'no continuation, \mu/10, restart',...
    'no continuation, \mu/10','no continuation, \mu',...
    'no continuation, 10 \mu','regular continuation','accelerated continuation'};
hl = legend(handles(cList+cOffset),names{cList+cOffset});
    
% legend(handles(cList+cOffset),names{cList+cOffset});

xlabel('iterations','fontsize',16);
ylabel('error','fontsize',16);
set(gca,'fontsize',16);
set(hl,'fontsize',20);
set(hl,'location','south')
ylim([1e-10,10]);
xlim([0,10001]);
%%
a = -.2;
orient landscape
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf,'PaperPosition',[a a (11-a) (8.5-a)]);
print('-dpdf','continuation_October_new2');
%%
% xlim([0,Inf])
%% add some info in the title
title(sprintf('Dantzig Problem, N=%d,M=%d,DCT matrix, nnz(x)=%d,x0=0\n||x||_\\infty=%.1f,\\mu=%.2e,\\mu for continuation = %.1f*\\mu,||D||=%.1e',...
    N,M,nnz(xExact),norm(xExact,Inf),mu0,mu_continuation/mu0,norm(D)) );
