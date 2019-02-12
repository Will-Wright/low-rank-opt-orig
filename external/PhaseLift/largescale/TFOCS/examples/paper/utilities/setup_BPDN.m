function [xNew,Anew,bNew,normA2,EPS,x0,lambda,D,omega] = setup_BPDN(opts)
% [xNew,Anew,bNew,normA2,EPS,xOriginal,lambda,D,omega] = setup_BPDN(opts)
%{
    Sets up A, b, ... for a basis pursuit denoising problem 

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
    TFOCS version 1.0, code by Michael Grant (mcg@cvxr.com) and Stephen Becker (srbecker@caltech.edu)

%}
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

SNR = setValue('SNR',30);

% either compressible or exactly sparse
SPARSE = setValue('sparse',true);

if SPARSE
    N = setValue('N',4*128);
    M = setValue('M',round(.5*N));  % want to have exact recovery when noiseless
    K = setValue('K',round(.25*M));
    
    DNR = setValue('DNR',20);
    x0        = zeros(N,1);
    omega     = sort(randsample(N,K));
    x0(omega) = sign(randn(K,1)).*10.^(DNR*rand(K,1)/20); % bug: used randn instead of rand until 6/25/10
else
    im = imread('cameraman.tif');
    waveletType = 'db1'; % same as Haar
    levels = 3;

    [c,s] = wavedec2(im,levels,waveletType);  % see 'wfilters' for options. Needs wavelet toolbox

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
    ORTH = setValue('orth',false); % make A orthogonal or not
    A   = randn(M,N);
    if ORTH
        ISOMETRY = true;
        A = orth(A')';
        pinvA = A';
        Apinv = @(b) A'*b;
        normAtA = 1;
    else
        ISOMETRY = false;  % i.e. pinv(A) ~= A'
        pinvA = pinv(A);
        Apinv = @(b) pinvA*b;
        normAtA = norm(A*A');
    end
    Af  = @(x) A*x;
    At  = @(y) A'*y;

elseif strcmpi( TYPE, 'DCT' )
    IMPLICIT = true;
    % Randomly pick rows of the sampling operator
    omega = randperm(N); omega = sort( omega(1:M) );
    downsample = @(x) x(omega,:);
    SS.type = '()'; SS.subs{1} = omega; SS.subs{2} = ':';
    upsample = @(x) subsasgn( zeros(N,size(x,2)),SS,x);
    
    Af = @(x) downsample(idct(x));
    At =@(x) dct(upsample(x));

    ISOMETRY = true;
    Apinv = At;
    normAtA = 1;
else
    error('bad value for ''type'' field in options');
end

fprintf('  USING N = %d, M = %d, MEASUREMENT MATRIX IS %s\n',...
    N, M, TYPE );

bExact    = Af(x0);
sigma     = 10^(-SNR/20)*(norm(bExact)/sqrt(M));
b_orig    = bExact + sigma * randn(M,1);
b         = b_orig;

% EPS       = sqrt( M + 2*sqrt(2*M) )*sigma;
EPS       = sqrt(M)*sigma;

% For sigma=0, often the l1 solution is the same as the original sparse solution
xExact      = x0;

%% Crude solution via 1st order algorithm
disp(' ---- Finding approximate solution to initial problem ---- ');
xPlug        = Apinv(b);           % initial guess
alpha = norm(xPlug,1);
beta = norm(xPlug)^2/2;
mu0 = (alpha/.95 - alpha )/beta;

if ~SPARSE
    mu0 = mu0/5;
end
mu = mu0;

xPlug     = zeros(N,1);          % initial guess
er        = @(x) norm(x-xExact)/norm(xExact);
er_plug   = @(x) norm(x-xPlug)/norm(xPlug);

oExact    = norm(xExact,1);
errob = @(f,l,xk) abs( norm(xk,1) - oExact ) / oExact; % xk not feasible so this may go negative
errp  = @(f,l,xk) norm( xk - xExact, 2 ) / norm( xExact,2 );
infeas = @(f,l,xk) (norm( Af(xk) - b )-EPS)/norm(b);

% Perform a preliminary solve
fprintf( '-----------\n' );
opts            = [];
opts.maxits     = 3000;
opts.printEvery = 500;
opts.tol        = 1e-1;
opts.errFcn     = {errob, errp, infeas};
% (OLD?) BUG: if L is provided but Lexact is NOT, then it fails!!! FIX.
% opts.L          = normAtA/mu;
opts.Lexact     = normAtA/mu;
opts.xPlug      = xPlug;
x               = xPlug;
opts.lambda0    = [];
opts.stopCrit = 2;  % legacy option

T = 1:N;
cnt = 0;
while (length(T) > M && cnt < 5) || (~EPS && (cnt < 4 && length(T) ~= K) )
    cnt = cnt + 1;
    xOld = x;
%     [x,lambda]=solve_BPDN( {Af,At}, b, mu, EPS, opts );
    if EPS
        [x,out]=solver_sBPDN( linop_handles([M,N],Af,At), b, EPS,mu,opts.xPlug,opts.lambda0,opts );
    else
        [x,out]=solver_sBP( linop_handles([M,N],Af,At), b,mu,opts.xPlug,opts.lambda0,opts );
    end
    lambda = out.dual;
    fprintf('Recovered solution is %.2e from original signal\n', norm(x-xExact)/norm(xExact) );
    fprintf('  and has %d nonzeros (N=%d, M=%d); original signal had %d nonzeros\n',...
        nnz(x),N,M,nnz(xExact) );
    T = find(x);
    opts.xPlug = x + cnt/(cnt+3)*( x - xOld );
    opts.lambda0 = lambda;
end

% solve once more, this time to high accuracy
% opts.tol = 1e-3;
opts.tol = 1e-8;
% [x,lambda]=solve_BPDN( {Af,At}, b, mu, EPS, opts );
if EPS
    [x,out]=solver_sBPDN( linop_handles([M,N],Af,At), b, EPS,mu,opts.xPlug,opts.lambda0,opts );
else
    [x,out]=solver_sBP( linop_handles([M,N],Af,At), b,mu,opts.xPlug,opts.lambda0,opts );
end
lambda = out.dual;
T = find(x);
fprintf('Recovered solution is %.2e from original signal\n', norm(x-xExact)/norm(xExact) );
fprintf('  and has %d nonzeros (N=%d, M=%d); original signal had %d nonzeros\n',...
    nnz(x),N,M,nnz(xExact) );


if length(T) > M
    disp('Warning: Failed to find vertex, trying to proceed after hard-thresholding');
%     error('Failed to find vertex after 10 continuation steps');
    xAbs = sort(abs(x),'descend');
    cutoff = xAbs(M);
    T = find( abs(x) > cutoff ); % should be < M strictly
end

%% Fix-up the crude solution to get a problem with known exact answer
% (i.e. x and y (y is aka lambda) are the exact solutions to a pertubed
% problem ).
% If EPS == 0, we'll modify A and x, and keep b and lambda fixed
% if EPS >  0, we'll modify A and b, and keep x and lambda fixed
Aty = At(lambda);
% KKT: need ||Aty||_inf <= 1, and on supp(x) it should be -sign(x)
% T = find(x); % see above now
Tc = setdiff(1:N,T);
% c = cond(A(:,T));  % remove this????
% if length(T) > M
%     error('PROBLEM: did not recover a vertex, cannot create exact solution');
% elseif c > 1e10
%     disp('PROBLEM: A is ill-conditioned on the support of x');
% end

xNew = x;

% == for EPS=0, we do things a little differently ==
if EPS == 0
    if IMPLICIT
        down = @(x) x(T,:);
        ST.type = '()'; ST.subs{1} = T; ST.subs{2} = ':';
        up = @(x) subsasgn( zeros(N,size(x,2)),ST,x);

        AA = implicitWrapper( @(x) Af(up(x)), @(y) down(At(y)) );
        [xNew(T),flag,relres,iter] = lsqr(AA,b, 1e-15, 500,[],[],x(T) );
        fprintf('  Cleaning up via LSQR; flag is %d, relres is %.2e, took %d iterations\n',flag,relres,iter);
    else
        xNew(T) = A(:,T)\b;
    end
    if ~all( sign(x) == sign(xNew) )
        fprintf('  Cleaned up primal version has different sign pattern though! Bad!\n');
    end
end
%%
% == Get a dual feasible variable by modifing A <-- A*D, D diagonal ==
D = ones(N,1);
r = .999; % require that ||(D*A'*y)_{T^c}||_\inf <= r < 1
D(Tc) = 1./max( abs(Aty(Tc))/r, 1);
% norm( D(Tc).* Aty(Tc), Inf )

D(T) = -1*sign(xNew(T))./Aty(T); % Sep, October: I think this is wrong. No, it's right.
% D(T) = 1./(sign(xNew(T).*Aty(T) )); % October: this fix isn't right: it gives all ones of course!
% D(T) = 1./abs(xNew(T).*Aty(T) ); % October: this fix isn't right: it gives all ones of course!
% (Note the negative: my lambda is like -lambda)

% if any( D < 0 ) 
%     error('  setup_BPDN: cleaning up solution for exact answer did not work!');
% end
if IMPLICIT
    Af = @(x) Af(x.*D); % if x is a matrix, this won't work (need repmat in that case, or make D diagonal)
    At = @(x) D.*At(x);

    ISOMETRY = false; % if diag(D) is not the identity, no longer have an isometry
    Apinv = [];  % not an easy calculation any more
%     normAtA = 1;
    [normA2,cnt] = my_normest( Af, At, N, 1e-6,200 ); % should be nearly 1
    normA2 = normA2^2;
    if cnt == 200
        disp(' warning: estimate of ||AA''|| not within desired tolerance');
    end
    
    % for output, give it in this form:
    Anew = {Af, At };
else
    Anew = A*diag(D);
    Af = @(x) Anew*x;
    At = @(x) Anew'*x;
    normA2 = norm(Anew*Anew');
end

if EPS == 0
% Need to re-solve for x, since now we need (A*D)*x = b
% This is trivial since D is diagonal. But we need D > 0,
% otherwise it means x changes sign, so we'd need to re-solve for y.
    xNew(T) = xNew(T)./D(T);
end

% == Get a primal feasible variable ==

if EPS
    z = -EPS*lambda/norm(lambda); % I have -lambda, so use the minus sign
    bNew = z + Af(xNew);
else
    bNew = b;
end


fprintf('  Modifying A <-- A.*diag(D). Entries of D: mean = 1+%.2e, var = %.2e\n',...
    mean(D)-1, var(D) );
if EPS
    fprintf('  Modifying b <-- b + z.  ||z|| = %.2e\n', norm( bNew - b )/norm(b) );
else
    fprintf('  Modifying x <-- xNew.  ||change|| = %.2e\n', norm( xNew - x )/norm(x) );
end

% the -1 in front of lambda because I have -lambda
pdgap = ( -bNew'*lambda - EPS*norm(lambda) ) - norm(xNew,1 );
fprintf('  Duality gap for our solutions is %.2e\n', pdgap );

oExactNew    = norm(xNew,1);

%% Verify via CVX
if N <= 1024 && ~IMPLICIT
    
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
        norm(bNew-Anew*x_cvx) <= EPS : lambda_cvx
cvx_end
% fprintf('Error is %.2e\n', norm(x_cvx - xExact)/norm(xExact) );
fprintf('CVX: rel. error in x is    %.2e\n', norm(x_cvx - xNew)/norm(xNew) );
f = @(x) norm(x,1);
fprintf('CVX: rel. error in f(x) is %.2e\n', (f(x_cvx)-f(xNew))/f(xNew));

if norm(xExact, Inf ) < 100*eps
    disp('-------------- CVX found all zero solution ---------');
end

end

%% for fun, see if we reconstructed the image well
if ~SPARSE && N == length(c) && false
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

