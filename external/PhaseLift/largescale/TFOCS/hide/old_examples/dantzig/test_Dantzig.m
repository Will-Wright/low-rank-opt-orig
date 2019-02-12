%{
    Run tests on the Dantzig Selector

Fri, June 25 2010

To run this file, use matlab's "cell" mode:
    Pick a "test" and execute that cell,
    then run the cell that calls "setup_BPDN" to generate that test.
    Each solver is in it's own cell, so call the appropriate cell
    To see the results, run the cell at the end of the file called "plot error_history"

%}

% modify this variable to point to the base directory
% of dropbox:
base = '~/Dropbox';

if ~exist('runOnce','var') || ~runOnce
    addpath( fullfile(base,'CNest') );
    addpath( fullfile(base,'CNest','examples','dantzig') );
    addpath( fullfile(base,'CNest','utilities') );
%     addpath( fullfile(base,'CNest','experiments','utilities') );
    %addpath ~/Documents/MATLAB/CNest/dantzig/
    %addpath ~/Documents/MATLAB/CNest/dantzig/utilities/
    runOnce = true;
end

%% Setup problem
%% Test 1
test=[];
test.name = 'test 1';
test.seed = 1232383;
test.SNR = 30; 
test.DNR = 20;
test.sparse = true;
N = 4*128;
M = round(.5*N);
K = round(.25*M );
test.N = N; test.M = M; test.K = K;
test1 = test;

%% Test 2  not k-sparse
test=[];
test.name = 'test 2';
test.seed = 83432;
test.SNR = 20; 
test.DNR = 60;
test.sparse = true;
N = 2^13;
M = round(.5*N);
% K = round(.25*M );
K = N;
test.N = N; test.M = M; test.K = K;
test.type = 'dct';
test1 = test;
%% Test 3
% run test 1, and then following modifications
test = test1;
test.name = 'test 3';
test.SNR = Inf;
% Used mu0 = (alpha/.95 - alpha )/beta; intermediate tol: 1e-6
%% Test 4
test=[];
test.name = 'test 4';
test.seed = 12321234;
test.SNR = 40; 
test.DNR = 10;  % a larger DNR does not mean exact soln has large DNR -- in fact, the opposite!
test.sparse = true;
N = 2^14;
M = round(N/8);
K = round(M/5 );
test.N = N; test.M = M; test.K = K;
test.type = 'dct';
%% Test 5  (Haar wavelets of cameraman image)
test=[];
test.name = 'test 5';
test.seed = 9234;
% test.SNR = Inf;   % hard to get a vertex solution with no noise
test.SNR = 30;
test.sparse = false;
% N = 2^13;  
% M = round(N/5);
% test.N = N; 
% test.M = M;
test.type = 'dct';
%% Test 6 Noiseless, N = 1024
test=[];
test.name = 'test 6';
test.seed = 93825;
test.SNR = Inf;
test.DNR = 10;
test.sparse = true;
N = 2^10;
M = round(.5*N);
K = round(.25*M );
test.N = N; test.M = M; test.K = K;
test.type = 'gaussian';
%{
Recovered solution is 1.61e-06 from original signal
  and has 128 nonzeros (N=1024, M=512); original signal had 128 nonzeros
  Modifying A <-- A.*diag(D). Entries of D: mean = 1+9.36e-07, var = 4.71e-10
  Modifying x <-- xNew.  ||change|| = 6.53e-05
  Duality gap for our solutions is 0.00e+00
CVX: rel. error in x is    4.54e-10
CVX: rel. error in f(x) is 4.03e-10  (SeDuMi)

CVX: rel. error in x is    3.70e-07
CVX: rel. error in f(x) is 3.31e-07 (SDPT3)

%}
%% Apply the test:
fprintf('Setting up exact solution for %s\n', test.name );
[xExact,A,b,normA2,EPS,D,xOrig,xPlug_exact,mu_exact] = setup_Dantzig(test);
delta0 = EPS;
if iscell(A)
    Aff = A{1};
    Att = A{2};
else
    Aff = @(x) A*x;
    Att = @(x) A'*x;
end
N = length(xExact); M = length(b);

Af = @(x) multiplyA( Aff, Att, 'forward', x );
At = @(x) multiplyA( Aff, Att, 'transpose', x );

% %% View dynamic range
% figure(2); clf;
% semilogy( sort(abs(xNew),'descend') );
% hold all
% semilogy( sort(abs(xOrig),'descend') );


%% Solve via conic solver

% if you want continuation, set it here:
% CONTINUATION = false;   % no continatuation
% CONTINUATION = 1;       % "standard" continuation
CONTINUATION = 2;         % "accelerated" continuation (recommended)

% xPlug = Att(b)/normA2;
% alpha = norm(xPlug,1);   beta = norm(xPlug)^2/2;
% want equal weights roughly, so try to solve for
% .5 = alpha / (alpha+mu*beta)
% mu0 = (alpha/.9 - alpha )/beta;
% mu0 = 1e-4;

if CONTINUATION == 2
%     mu0 = .1*mu0;  %10000 for test 1
    mu0 = 1e-3;
elseif CONTINUATION == 1
    mu0 = 1;
end
% mu0 = 1e-2;
% mu0 = 1e-1;

mu=mu0;


% er        = @(x) norm(x-xExact)/norm(xExact);
% er_plug   = @(x) norm(x-xExact)/norm(xExact);

% In these error computations, do NOT call multiplyA since we do NOT want
% to count these.
oExact    = norm(xExact,1);
errp  = @(f,l,xk) norm( xk - xExact, 2 ) / norm( xExact,2 );
% infeas = @(f,l,xk) (norm( D.*Att(Aff(xk) - b) )-EPS)/norm(b);

globalErrors = @(f,l,xk) ComputeErrors( xExact, oExact, xk, f );
multiplyA();  % zero it out.

opts = [];
opts.errFcn     = {errp, globalErrors};
opts.Lexact     = normA2/mu;
x = zeros(N,1);
% x = xPlug_exact;  % !!!
xOld = x;

% tol = 1e-2*norm(bExact);
tol = 1e1*norm(b);
opts.tol = tol;
if CONTINUATION
    k_end = 14;
else
    k_end = 1;  % no continuation
end

for k = 1:k_end
    if CONTINUATION == 2
        opts.x0 = x + (k/(k+3))*(x-xOld); % accelerated continuation
        xOld = x;
    else
        opts.x0 = x;
    end
    % at the final outer iteration, ask for more precision:
    if k == k_end
        disp('  == FINAL ITERATION == ');
%         mu =  mu0/10;
        opts.tol = tol*1e-6;
%         opts.Lexact     = normA2/mu;
    else
        opts.tol = max( opts.tol/2, tol*1e-8 );
    end
    [x,lambda,out]=solve_Dantzig( {Af,At}, b,D, mu, EPS, opts );
    opts.lambda0 = lambda;
end

error_history = multiplyA();  % record error history

%% plot error_history
% May not be smooth, due to backtracking, etc.
figure(1);
iter = 1:length(error_history);
semilogy( iter, error_history(:,2 ) );
hold all
xlim([0,3000])
% xlim([0,Inf])
%{
    Errors:
column 1: iteration count
column 2: norm(x-xTrue)/norm(xTrue);
column 3: norm(x-xTrue,'inf');
column 4: length(T1)+length(T2)-2*length(intersect(T1,T2) );
column 5: abs( norm(x,1) -fTrue)/abs(fTrue);
%}
%% view solution
% figure(2);
% hold all;
% semilogy( sort(abs(x),'descend') );
