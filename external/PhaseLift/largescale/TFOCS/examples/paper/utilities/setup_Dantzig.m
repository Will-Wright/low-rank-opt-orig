function [xNew,Anew,bNew,normA2,delta0,D,x0,xPlug,mu,omega,lambdaExact] = setup_Dantzig(opts)
% [xNew,Anew,bNew,normA2,delta0,DNew,xOriginal,xPlug,mu,omega,lambda] = setup_Dantzig(opts)
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


see "perturbToExact", a better replacement for getting exact solution

TFOCS ver

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

setValue = @(field,default)setValue3arg(field,default,opts);

% Random number generator
MersenneTwister = 'mt19937ar'; 
seed = setValue('seed', sum(10000*clock) );
if verLessThan('matlab','7.7')
    randn('state',seed);
    rand('state',seed);
else
    RandStream.setDefaultStream(RandStream(MersenneTwister,'seed',seed) );
end

SMOOTHED = setValue( 'smoothed', false );
SNR = setValue('SNR',30);

% ----------------- Setup the original signal ----------------------------
% either compressible or exactly sparse
SPARSE = setValue('sparse',true);
if SPARSE
    % The original signal is k-sparse
    N = setValue('N',512);
    M = setValue('M',round(.5*N));  % want to have exact recovery when noiseless
    K = setValue('K',round(.25*M));
    
    DNR = setValue('DNR',20);
    x0        = zeros(N,1);
    omega     = randperm(N);
    omega     = omega(1:K);
    x0(omega) = sign(randn(K,1)).*10.^(DNR*rand(K,1)/20);
else
    % The original signal is the wavelet transform of an image
    im = imread('cameraman.tif'); % this image comes standard with MATLAZB
    waveletType = 'db1'; % aka Haar Wavelet
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

% ----------------- Setup the sampling matrix ----------------------------
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
        At =@(x) idct(upsample(x));
        Apinv = At;
    end
else
    error('bad value for ''type'' field in options');
end

fprintf('  USING N = %d, M = %d, MEASUREMENT MATRIX IS %s',...
    N, M, upper(TYPE) );
if SPARSE, fprintf(', sparsity is %d, dynamic range is %.1f dB', K, DNR ); end
fprintf('\n');

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


% -------- estimate ||DA'A||^2 --------------------------------------------
if IMPLICIT
    temp = linop_handles( [N,N], @(x)D.*(At(Af(x))), @(y)At(Af(D.*y)), 'R2R' );
    [normA2,cnt] = linop_normest( temp, 'R2R', 1e-5, 800 );
    normA2 = normA2^2;
    if cnt == 600
        disp(' warning: estimate of ||AA''|| not within desired tolerance');
        normA2 = (1.03)*normA2;  % be on the safe side
    end
else
    normA2 = norm( diag(D)*(A'*A) )^2;
end



% For sigma=0, often the l1 solution is the same as the original sparse solution
xExact      = x0;

% -------- "Crude" solution via 1st order algorithm -----------------------
disp(' ---- Finding approximate solution to initial problem ---- ');
xPlug        = Apinv(b);           % initial guess
% mu        = .00001*norm(x0,'inf'); % smoothing parameter
alpha = norm(xPlug,1);
beta = norm(xPlug)^2/2;
mu0 = (alpha/.95 - alpha )/beta;

if ~SPARSE
    mu0 = mu0/5;
end
mu0 = 20*mu0;

if isfield(opts,'mu') && ~isempty(opts.mu) && opts.mu > 0
    FORCE_MU = true;
    MU = opts.mu;
else
    FORCE_MU = false;
end

mu = mu0;

% xPlug     = zeros(N,1);          % initial guess
xPlug     = setValue('xPlug',zeros(N,1)); 
er        = @(x) norm(x-xExact)/norm(xExact);
er_plug   = @(x) norm(x-xPlug)/norm(xPlug);

oExact    = norm(xExact,1);
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
% opts.lambda0    = [];
lambda          = [];
opts.maxits     = 3000;
opts.printEvery = 500;
tol = 1e-2*norm(b);
opts.tol = tol;
opts.errFcn     = {errob, errp, infeas};
opts.xPlug      = xPlug;
x               = xPlug;
xOld            = x;

% NewFormat       = true;  % apply changes made on Sep 8 or not
NewFormat = setValue('NewFormat',false);
fprintf('Using new format for perturbing to exact solution\n');

T = 1:N;
cnt = 0;
if N > 2^10
    cnt_max = 20;
else
    cnt_max = 12;
end
high_tol = 1; % added Sep 8 2010
if SMOOTHED
    cnt_max = 0; 
    mu = 10*mu;
    if NewFormat, high_tol = 1e4; end
end  % do NOT do continuation
cnt_max=setValue('cont',cnt_max);   % override if value is given

mu0 = mu;
while cnt <= cnt_max
    mu = mu0;
    opts.tol = max( opts.tol/2, tol*1e-8 );
    if cnt == cnt_max 
        opts.tol    = tol*1e-8;
        opts.tol    = opts.tol / high_tol; % added Sep 8
        mu          = mu0/100;
        opts.maxits = 4000; % Sept 8, trying below:
        if SMOOTHED && NewFormat
            opts.maxits = 12000; % also Sep 8
            opts.restart = 1000; % helps!
        end
    end
    if FORCE_MU, mu = MU; end  % override mu if the option has been set
    opts.Lexact = normA2/mu;
    cnt = cnt + 1;
    opts.xPlug = x + cnt/(cnt+3)*( x - xOld );
    xOld = x;
%     [x,lambda]=solve_Dantzig( {Af,At}, b, D, mu,delta0, opts );
    A = linop_handles([N,N],Af,At);
    opts.stopCrit = 2;  % legacy option
    [ x, out ] = solver_sDantzig( {A,D}, b, delta0, mu,opts.xPlug,lambda, opts );
    lambda = out.dual;

    T = find(x);
    opts.lambda0 = lambda;
    if cnt <= cnt_max, disp('   ---------- continuation --------- ' ); end
end
fprintf('Recovered solution is %.2e from original signal\n', norm(x-xExact)/norm(xExact) );
fprintf('  and has %d nonzeros (N=%d, M=%d); original signal had %d nonzeros. mu is %.2e\n',...
    nnz(x),N,M,nnz(xExact), mu );


% for the non-smoothed case, we want a LP, so solution should be at vertex
if length(T) > M && ~SMOOTHED 
    disp('Warning: Failed to find vertex, so hard-thresholding');
    xAbs = sort(abs(x),'descend');
    cutoff = xAbs(M);
    T = find( abs(x) >= cutoff ); % should be < M strictly
    
    % Sep 8 2010, adding this:
    x( setdiff(1:N,T) ) = 0;
end
x_tfocs = x;


% ---- Fix-up the crude solution to get a problem with known exact answer -
% We will keep the dual solution lambda (or perhaps hard-threshold it),
% and keep A, xPlug, delta0 and mu.
% We perturb D and b, and generate a new primal solution.

% see the newer function "perturbToExact.m"

if numel(D) == 1, D = D*ones(N,1); end
xPlug       = opts.xPlug;
oldD        = D;
if SMOOTHED && NewFormat
    [xNew,bNew,D,lambdaExact,z] = perturbToExact(lambda,x,Af,At,D,delta0,b,mu,xPlug);
else

    if SMOOTHED && NewFormat
        lambdaExact = lambda .* ( abs(lambda) > 1e-5 * norm(lambda,Inf) );
    elseif NewFormat
        lambdaExact = lambda .* ( abs(lambda) > 1e-5 * norm(lambda,Inf) );
    else
        lambdaExact = lambda .* ( abs(lambda) > 1e-8 * norm(lambda,Inf) );
    end
    if delta0
        z           = b - Af(x);
        Atz         = At(z);    % added Oct 25
        DAtz        = D.*Atz;
        fprintf('  Initial x from solver is infeasible by %.2e\n',norm(DAtz,Inf) - delta0 );
%         Dexact      = D .* -sign(lambdaExact).*( delta0 ./ DAtz );
        Dexact      = sign(lambdaExact).*( delta0 ./ Atz ); % Oct 25, using this
        SGN         = -1; % Oct 25
        
        if SMOOTHED && NewFormat
            tt      = lambdaExact ~= 0;
        elseif NewFormat
            tt      = lambdaExact ~= 0;
        else
            tt      = lambda ~= 0;  % don't want to break continuation experiment
        end
        D(tt)       = Dexact(tt);
        tt          = ~tt & ( abs(DAtz) > delta0 );
        D(tt)       = 0.999 * Dexact(tt);
    else
        % delta0 = 0, and D is 1 (or anything really).  Need z = 0.
        z           = zeros(size(b));
        % do not modify D
    end
    
    if ~SMOOTHED && NewFormat
%         mu  = mu/100000;
        % get the x0 out of the picture
        xNew        = shrink( xPlug - SGN*At(Af(D.*lambdaExact))/mu, 1/mu );
    else
        xNew        = shrink( xPlug - SGN*At(Af(D.*lambdaExact))/mu, 1/mu );
    end
    bNew        = z + Af( xNew );
end

% Sep 8
if any( D < 0 )
    disp(' - - - - - warning: D < 0 - - - - - - - ');
elseif any( D == 0 )
    disp(' - - - - - warning: D_i = 0 for some i - - - - ');
end

% ------------ re-estimate ||DA'A||^2, since D has changed ---------------
if IMPLICIT
    temp = linop_handles([N,N], @(x)D.*(At(Af(x))), @(y)At(Af(D.*y)), 'R2R' );
    [normA2,cnt] = linop_normest( temp, 'R2R', 1e-6, 400 );
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
fprintf('  Modifying D.  ||D-D_old||_inf/||D_old||_inf is %.2e\n',...
    norm( D - oldD,Inf ) / norm(oldD,Inf ) );
fprintf('  Modifying b.  ||b-b_old||/||b_old|| is %.2e\n',...
    norm( bNew - b )/norm(b) );
fprintf('  Modifying x.  ||x-x_old||/||x_old|| is %.2e\n',...
    norm( xNew - x )/norm(x) );

pdgap1 = ( SGN*lambdaExact'*(D.*At(z)) + delta0*norm(lambdaExact,1) );
pdgap2 = ( -lambdaExact'*(D.*At(z)) + delta0*norm(lambdaExact,1) ); % Oct 25
% pdgap = min( abs(pdgap1), abs(pdgap2) ); % depends on our convention for lambda
pdgap = pdgap1;

fprintf('  Primal feasibility violated by: %.2e\n', norm( D.*At(z), 'inf') - delta0 );
fprintf('  Duality gap is %.2e;   mu is %.2e\n', pdgap,mu );


end % end of main function

% Helper function
function q=setValue3arg(field,default,opts)
 if ~isfield(opts,field)
    opts.(field) = default;
 end
 q = opts.(field);
end

