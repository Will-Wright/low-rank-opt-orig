function [xNew,Anew,bNew,normA2,delta0,D,x0,xPlug,mu] = setup_Dantzig(opts)
% [xNew,Anew,bNew,normA2,delta0,DNew,xOriginal] = setup_Dantzig(opts)
%    to setup a problem, or
% setup_Dantzig() to see available options
%{
    Sets up A, b, ... for a Dantzig Selector problem 

    You can control:
        M, N, K
        Dynamic Range (DNR)
        Signal-to-noise ratio (SNR)
        k-Sparse signal, or realistic signal (wavelet coefficients of an image)
        Random seed
        Type of measurement matrix (e.g. ~N(0,1) entries, or partial DCT)
        Format of measurement matrix: explicit matrix, or use a function handle

    This script will solve, then use that solution as the exact solution
    to a slightly perturbed problem.

    Stephen Becker, June 2010

    To see the available options, call this function without any output
    arguments.

    Might depend on hadamard mex file!

%}

if nargout == 0
    % Display available options
    disp('  setup_Dantzig(opts): possible fields for the "opts" structre: ');
    disp('      seed (scalar): seeds pseudo-random number generator')
    disp('      SNR (scalar, or Inf): controls Signal-To-Noise ratio');
    disp('      sparse (logical): make signal with K nonzeros, or (if false) use wavelet coefficients');
    disp('      M, N, K (integers): [M,N] = size(A), K = nnz(x)');
    disp('      DNR (scalar): dynamic range of x in dB');
    disp('      type (string): ''gaussian'' or ''dct'', for measurement matrix');
    disp('          or ''hadmard'' or ''fwht'' for Hadamard transform ');
    disp('      smoothed (logical): if false (default), sets up solution to unsmoothed problem');
    disp(' ');
    return;
end

if nargin == 0, opts = []; end

% for some reason, can't use nested functions when you need to call cvx :
% it gives a weird "Attempt to add "cvx_problem" to a static workspace"
setValue = @(field,default)setValue3arg(field,default,opts);

% Random number generator
MersenneTwister = 'mt19937ar'; % they make this cumbersome!
seed = setValue('seed', sum(10000*clock) );
if verLessThan('matlab','7.7')
    randn('state',seed);
    rand('state',seed);
else
    RandStream.setDefaultStream(RandStream(MersenneTwister,'seed',seed) );
end

SMOOTHED = setValue( 'smoothed', false );
SNR = setValue('SNR',30);

% either compressible or exactly sparse
SPARSE = setValue('sparse',true);

if SPARSE
    N = setValue('N',512);
    M = setValue('M',round(.5*N));  % want to have exact recovery when noiseless
    K = setValue('K',round(.25*M));
    
    DNR = setValue('DNR',20);
    x0        = zeros(N,1);
    omega     = randperm(N);
    omega     = omega(1:K);
    x0(omega) = sign(randn(K,1)).*10.^(DNR*rand(K,1)/20);
else
    im = imread('cameraman.tif');
    waveletType = 'db1'; % same as Haar
    levels = 3;

    [c,s] = wavedec2(im,levels,waveletType);  % see 'wfilters' for options. Needs wavelet toolbox
%     d = detcoef2( 'h', c, s, 1 );
%     d = [d,detcoef2( 'd', c, s, 1 )];
%     a = appcoef2( c, s, 'db3', 1);
%     d = [d;detcoef2( 'v', c, s, 1 ), a];
%     imshow(d);
%     semilogy( sort(abs(c),'descend') );

    N = setValue('N',length(c) );
    if N > length(c)
        error('Bad value for N when using wavelet coefficients');
    elseif N < length(c)
        x0 = c( randsample(length(c), N) ) ; % requires stats toolbox
    else
        x0 = c;
    end
    % To help make the wavelet coefficients less coherent with the DCT
    % measurement, we will randomly permute them:
    prm = randperm(N);
    x0 = x0(prm);
    
    x0 = x0';  % change row vector to column vector
    M = setValue('M',round(.3*N));

    K = N;
end

TYPE = setValue('type','GAUSSIAN');
if strcmpi( TYPE, 'GAUSSIAN' )
    IMPLICIT = false;
    A   = randn(M,N);
    Af  = @(x) A*x;
    At  = @(y) A'*y;
    pinvA = pinv(A);
    Apinv = @(b) pinvA*b;
elseif strcmpi( TYPE, 'DCT' ) || strcmpi(TYPE,'hadamard') ||strcmpi(TYPE,'fwht')
    IMPLICIT = true;
    % Randomly pick rows of the sampling operator
    omega = randperm(N); omega = sort( omega(1:M) );
    downsample = @(x) x(omega,:);
    SS.type = '()'; SS.subs{1} = omega; SS.subs{2} = ':';
    upsample = @(x) subsasgn( zeros(N,size(x,2)),SS,x);
    if strcmpi(TYPE,'hadamard') ||strcmpi(TYPE,'fwht')
        % can use MATLAB's "fwht" code, but it's slow
        % Much better is Peter's hadamard mex file
        
        % Add in a random permutation to get better incoherence:
        rp = randperm(N);
        [ignore,rp_inv] = sort(rp);
        rp_inv_f = @(x) x(rp_inv,:);
        
        Af = @(x) downsample(hadamard(x(rp,:)));
        At =@(x) rp_inv_f(hadamard(upsample(x)));
        Apinv = @(x)At(x)/N;
    else
        Af = @(x) downsample(dct(x));
        At =@(x) idct(upsample(x)); % 7/1/10 switching dct and idct
        Apinv = At;
    end
else
    error('bad value for ''type'' field in options');
end

fprintf('  USING N = %d, M = %d, MEASUREMENT MATRIX IS %s\n',...
    N, M, TYPE );

bExact    = Af(x0);
sigma     = 10^(-SNR/20)*(norm(bExact)/sqrt(M));
b_orig    = bExact + sigma * randn(M,1);
b         = b_orig;

% Construct the error bounds
alpha   = .05; % Problem will be feasible 100*(1-alpha) percent of the time
if IMPLICIT
%     Anorms = findColumnNorms(Af,N);
    Anorms = findColumnNorms(At,N,M);  % faster if M < N
    Anorms = Anorms';
else
    Anorms  = sqrt(sum(A.^2))';
end
nTrials = min(4*N,400); 
z       = randn(M,nTrials);
supAtz  = sort(max(At(z) ./Anorms(:,ones(1,nTrials))));
thresh  = supAtz(round(nTrials*(1-alpha)));  % empirical
delta   = thresh*sigma*Anorms;  % a vector, for Dantzig solvers
if all(delta > 1e-10)
    D       = 1./delta; % watch out for division by 0
    delta0  = 1;  % instead of using delta, we use D and delta0 = 1
else
    D       = 1;
    delta0  = 0;
end
% normA   = norm(A);

if IMPLICIT
    [normA2,cnt] = my_normest( @(x)D.*(At(Af(x))),...
        @(y)At(Af(D.*y)), N, 1e-6,400 ); % should be nearly 1
    normA2 = normA2^2;
    if cnt == 400
        disp(' warning: estimate of ||AA''|| not within desired tolerance');
        normA2 = (1.05)*normA2;  % be on the safe side
    end
else
    normA2 = norm( diag(D)*(A'*A) )^2;
end



% For sigma=0, often the l1 solution is the same as the original sparse solution
xExact      = x0;

%% Crude solution via 1st order algorithm
disp(' ---- Finding approximate solution to initial problem ---- ');
xPlug        = Apinv(b);           % initial guess
% mu        = .00001*norm(x0,'inf'); % smoothing parameter
alpha = norm(xPlug,1);
beta = norm(xPlug)^2/2;
% beta = norm(xPlug-yk)^2/2;
% want equal weights roughly, so try to solve for
% .5 = alpha / (alpha+mu*beta)
mu0 = (alpha/.95 - alpha )/beta;

if ~SPARSE
    mu0 = mu0/5;
end
% if delta0
%     mu0 = mu0/4;
% end
mu0 = 20*mu0;
mu = mu0;

xPlug     = zeros(N,1);          % initial guess
er        = @(x) norm(x-xExact)/norm(xExact);
er_plug   = @(x) norm(x-xPlug)/norm(xPlug);

oExact    = norm(xExact,1);
% errob = @(f,l,xk) abs( f - oExact ) / oExact;
errob = @(f,l,xk) abs( norm(xk,1) - oExact ) / oExact; % xk not feasible so this may go negative
errp  = @(f,l,xk) norm( xk - xExact, 2 ) / norm( xExact,2 );
if delta0 > 0
    infeas = @(f,l,xk) (norm( D.*(At(b-Af(xk))), 'inf' ) - delta0 )/delta0;
else
    infeas = @(f,l,xk) norm( D.*(At(b-Af(xk))), 'inf' );
end

% Perform a preliminary solve
fprintf( '-----------\n' );
opts            = [];
opts.maxits     = 3000;
opts.printEvery = 500;
tol = 1e-2*norm(b);
opts.tol = tol;
opts.errFcn     = {errob, errp, infeas};
% opts.L          = normA2/mu;
opts.Lexact     = normA2/mu;
opts.xPlug      = xPlug;
x               = xPlug;
xOld            = x;

T = 1:N;
cnt = 0;
if N > 2^10
    cnt_max = 20;
else
    cnt_max = 10;
end
if SMOOTHED
    cnt_max = 0; 
    mu = 10*mu;
end  % do NOT do continuation

while cnt <= cnt_max
% while (length(T) > M && cnt < 5) || (cnt < 4 && length(T) ~= K)
    opts.tol = max( opts.tol/2, tol*1e-8 );
    if cnt == cnt_max 
        opts.tol    = tol*1e-8;
        mu          = mu/100;
        opts.Lexact = normA2/mu;
        opts.maxits = 4000;
    end
    cnt = cnt + 1;
    opts.xPlug = x + cnt/(cnt+3)*( x - xOld );
    xOld = x;
    [x,lambda]=solve_Dantzig( {Af,At}, b, D, mu,delta0, opts );
    
    fprintf('Recovered solution is %.2e from original signal\n', norm(x-xExact)/norm(xExact) );
    fprintf('  and has %d nonzeros (N=%d, M=%d); original signal had %d nonzeros\n',...
        nnz(x),N,M,nnz(xExact) );
    T = find(x);
    opts.lambda0 = lambda;
end

% T = find(x);
% fprintf('Recovered solution is %.2e from original signal\n', norm(x-xExact)/norm(xExact) );
% fprintf('  and has %d nonzeros (N=%d, M=%d); original signal had %d nonzeros\n',...
%     nnz(x),N,M,nnz(xExact) );


if length(T) > M
    disp('Warning: Failed to find vertex, trying to proceed after hard-thresholding');
%     error('Failed to find vertex after 10 continuation steps');
    xAbs = sort(abs(x),'descend');
    cutoff = xAbs(M);
    T = find( abs(x) > cutoff ); % should be < M strictly
end

%% Fix-up the crude solution to get a problem with known exact answer
% We will keep the dual solution lambda (or perhaps hard-threshold it),
% and keep A, xPlug, delta0 and mu.
% We perturb D and b, and generate a new primal solution.
if numel(D) == 1, D = D*ones(N,1); end
xPlug       = opts.xPlug;
lambdaExact = lambda .* ( abs(lambda) > 1e-8 * norm(lambda,Inf) );
if delta0
    z           = b - Af(x);
    DAtz        = D.*At( z );
    oldD        = D;
    Dexact      = D .* -sign(lambdaExact).*( delta0 ./ DAtz );
    tt          = lambda ~= 0;
    D(tt)       = Dexact(tt);
    tt          = ~tt & ( abs(DAtz) > delta0 );
    D(tt)       = 0.999 * Dexact(tt);
else
    % delta0 = 0, and D is 1 (or anything really).  Need z = 0.
    z           = zeros(size(b));
    oldD        = D; % do not modify D
end
xNew        = shrink( xPlug - At(Af(D.*lambdaExact))/mu, 1/mu );
% bNew        = b + Af( xNew - x );
bNew        = z + Af( xNew );

obj = @(x) norm(x,1) + mu/2*norm(x-xPlug)^2;
oExact      = obj(xNew);

if IMPLICIT
    [normA2,cnt] = my_normest( @(x)D.*(At(Af(x))),...
        @(y)At(Af(D.*y)), N, 1e-6,400 ); % should be nearly 1
    normA2 = normA2^2;
    if cnt == 400
        disp(' warning: estimate of ||AA''|| not within desired tolerance');
        normA2 = (1.05)*normA2;  % be on the safe side
    end
    Anew = {Af,At};
else
    normA2 = norm( diag(D)*(A'*A) )^2;
    Anew = A;
end


fprintf('  Modifying D.  ||D-D_old||/||D_old|| is %.2e\n',...
    norm( D - oldD ) / norm(oldD ) );
fprintf('  Modifying b.  ||b-b_old||/||b_old|| is %.2e\n',...
    norm( bNew - b )/norm(b) );
fprintf('  Modifying x.  ||x-x_old||/||x_old|| is %.2e\n',...
    norm( xNew - x )/norm(x) );

pdgap = ( lambdaExact'*(D.*At(z)) + delta0*norm(lambdaExact,1) );

fprintf('  Primal feasibility violated by: %.2e\n', norm( D.*At(z), 'inf') - delta0 );
fprintf('  Duality gap is %.2e;   mu is %.2e\n', pdgap,mu );


%% Verify via CVX (for un-smoothed version)
if N <= 1024 && ~IMPLICIT
if delta0    
    cvx_begin
    % SeDuMi is default solver
    %     cvx_solver sdpt3
          cvx_precision best
    %     cvx_precision high
           cvx_quiet true
           variable x_cvx(N,1)
           dual variable lambda_cvx
           minimize norm( x_cvx, 1 )
           subject to
                abs(diag(D)*A'*(bNew-A*x_cvx)) <= delta0 : lambda_cvx
    cvx_end
else 
    cvx_begin
    % SeDuMi is default solver
    %     cvx_solver sdpt3
          cvx_precision best
    %     cvx_precision high
           cvx_quiet true
           variable x_cvx(N,1)
           dual variable lambda_cvx
           minimize norm( x_cvx, 1 )
           subject to
                bNew == A*x_cvx : lambda_cvx
    cvx_end
end




% fprintf('Error is %.2e\n', norm(x_cvx - xExact)/norm(xExact) );
fprintf('CVX: rel. error in x is    %.2e\n', norm(x_cvx - xNew)/norm(xNew) );
f = @(x) norm(x,1);
fprintf('CVX: rel. error in f(x) is %.2e\n', (f(x_cvx)-f(xNew))/f(xNew));

if norm(xExact, Inf ) < 100*eps
    disp('-------------- CVX found all zero solution ---------');
end


% -- now, have CVX solve smoothed problem --
if delta0    
    cvx_begin
    % SeDuMi is default solver
    %     cvx_solver sdpt3
          cvx_precision best
    %     cvx_precision high
           cvx_quiet true
           variable x_cvx(N,1)
           dual variable lambda_cvx
           minimize norm( x_cvx, 1 ) + mu/2*sum_square(x_cvx-xPlug)
           subject to
                abs(diag(D)*A'*(bNew-A*x_cvx)) <= delta0 : lambda_cvx
    cvx_end
else 
    cvx_begin
    % SeDuMi is default solver
    %     cvx_solver sdpt3
          cvx_precision best
    %     cvx_precision high
           cvx_quiet true
           variable x_cvx(N,1)
           dual variable lambda_cvx
           minimize norm( x_cvx, 1 ) + mu/2*sum_square(x_cvx-xPlug)
           subject to
                bNew == A*x_cvx : lambda_cvx
    cvx_end
end


% fprintf('Error is %.2e\n', norm(x_cvx - xExact)/norm(xExact) );
fprintf('CVX on smoothed problem: rel. error in x is    %.2e\n', norm(x_cvx - xNew)/norm(xNew) );
f = @(x) norm(x,1);
fprintf('CVX on smoothed problem: rel. error in f(x) is %.2e\n', (f(x_cvx)-f(xNew))/f(xNew));

if norm(xExact, Inf ) < 100*eps
    disp('-------------- CVX found all zero solution for smoothed problem ---------');
end

end

%% for fun, see if we reconstructed the image well
if ~SPARSE && N == length(c)
    figure(2); clf;
%     imshow( im );

    % undo the permutation
    [ignore,prmReverse] = sort(prm);
    xx = x(prmReverse);
    im_hat = waverec2(xx',s,waveletType);
    imshow( [im, im_hat] );
    
end

end % end of main function

function q=setValue3arg(field,default,opts)
  if ~isfield(opts,field)
      opts.(field) = default;
  end
  q = opts.(field);
end

