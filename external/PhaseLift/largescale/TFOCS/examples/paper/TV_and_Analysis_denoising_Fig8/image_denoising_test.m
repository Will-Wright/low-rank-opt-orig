%{
    Example of Total Variation norm (TV) recovery
    combined with analysis (e.g. weighted l1)

    Recreates the figures in the paper

    Should take minute or so to run this file.

    Stephen Becker, July 2010   srbecker@caltech.edu    

Toolboxes required:
communication_toolbox (awgn)
image_toolbox  (imshow, dct)
wavelet_toolbox (for the wavelets) -- if you don't have this,
    code will use dct, which requires image_toolbox
%}

% Try to make some simple replacement functions
% for users who don't have the toolboxes

% clear wavedec2 waverec2 imshow awgn
if ~license('test','communication_toolbox')
    awgn = @(x,snr,ignore) x + ...
        10^( (10*log10(sum(abs(x(:)).^2)/length(x(:))) - snr)/20 )*randn(size(x));
end
% If you don't have the image toolbox, must define imshow and the 2D DCT
if ~license('test','image_toolbox')
    imshow  = @my_imshow;
    dct2    = @(x)dct(dct(x).').';
end

% If you don't have the wavelet toolbox, use the DCT instead
if ~license('test','wavelet_toolbox')
    wavedec2 = @(varargin) vec(dct2(mat(varargin{1})));
    waverec2 = @(varargin) vec(idct2(mat(varargin{1})));
end


%%

%%%%%%%%%%%%%%%
% Clean image %
%%%%%%%%%%%%%%%

X_exact = double(imread('cameraman.tif'))/255;
szX     = size(X_exact);
N       = prod(szX);
mat     = @(X) reshape(X,szX);
vec     = @(X) reshape(X,N,1);
maxI    = 1;
PSNR    = @(X) 20*log10(maxI*sqrt(N)/norm(mat(X)-X_exact,'fro')); 
TV      = @(X) norm([diff(mat(X),1,2),zeros(n1,1)]+1j*[diff(mat(X),1,1);zeros(1,n2)],'fro');

%%%%%%%%%%%%%%%
% Noisy image %
%%%%%%%%%%%%%%%

randn('state',sum(double('denoising test')) );
SNR     = 20;
X_noisy = awgn( X_exact, SNR, 'measured' );
noise   = norm( X_noisy - X_exact, 'fro' ) / sqrt(N);
fprintf( 'Adding noise of std %.2f for a SNR of %.1f dB\n', noise, SNR );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wavelet thresholding oracle %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Periodization extension mode
% Symmetric 9/7 biothogonal wavelets---what JPEG-2000 uses
% See 'wfilters' for options

waveletType     = 'bior4.4';
transposeType   = ['rbio', waveletType(end-2:end) ];

if ~license('test','wavelet_toolbox')
    nLevels = 1; % not used
    c = wavedec2( X_noisy, nLevels, waveletType );
    s = [];
else
    dwtmode('per');
    nLevels = min( 3, wmaxlev(szX,waveletType) );
    [c,s] = wavedec2( X_noisy, nLevels, waveletType );
end

% Try a range of thresholds to find the best PSNR
cSort = sort( abs(c), 'descend' );
errBest = -Inf;
ind_levels = cSort .^ 2;
ind_levels = cumsum(ind_levels) ./ sum(ind_levels);
for gamma = .98:.001:.999
    ind = find( ind_levels > gamma, 1 );
    if isempty(ind), ind = N; end
    c2     = c .* ( abs(c) > cSort(ind) );
    X_test = waverec2(c2,s,waveletType);
    err    = PSNR(X_test);
    if err > errBest, 
        errBest = err;
        gammaBest = gamma;
        X_oracle = X_test;
    end
end
X_oracle = mat(X_oracle);

%%%%%%%%%%%%%%%%%%%%%%%
% Show results so far %
%%%%%%%%%%%%%%%%%%%%%%%

imshow( [X_exact,X_noisy,X_oracle] );
title( sprintf('Original\t\tNoisy (PSNR %.1f dB)\tOracle Wavelet Thresholding (PNSR %.1f dB)',...
    PSNR(X_noisy), PSNR(X_oracle) ) );

%%%%%%%%%%%%%%%%%%%%
% Wavelet analysis %
%%%%%%%%%%%%%%%%%%%%

Wf  = @(X) vec(wavedec2(mat(X),nLevels,waveletType));
Wt  = @(c) vec(waverec2(vec(c),s,transposeType));
W_wavelet       = linop_handles([N,N], Wf, Wt );
normWavelet     = linop_test(W_wavelet);
%%

% Measurements:
M = N;
n1 = szX(1);
n2 = szX(2);
A = linop_handles({[n1*n2,1],[M,1]}, @(x)x, @(x)x );
EPS = norm( X_noisy - X_exact, 'fro' );
Xbp = X_noisy; % backprojection solutution, since pinv(A) = I
xbp = vec(Xbp);
b   = xbp;
W_tv    = linop_TV( szX );
normTV           = linop_TV( szX, [], 'norm' );
opts             = [];
opts.stopCrit = 2;  % legacy option
opts.normW12     = normTV^2;
opts.normW22     = normWavelet^2;
opts.normA2      = 1; % "A" is the identity

X_recovered = {};
%% Set some universal parameters
% k_end is the max number of continuation steps.
% For this test, it's not a big deal whether we use continuation or not
% (since we are denoising, we have a pretty good guess for x0,
%  so the ||x - x_0 ||^2 term is not very significant )
k_end = [1,1,1];

% max # if iterations
maxits = [100,100,100];
%% for both tv and wavelets
testNo = 1;
string{testNo} = 'TV + Analysis';
% beta   = 5; alpha  = 1;  mu = 160; % value in paper. 
beta   = 1; alpha  = 1;  mu = 25; % visually pleasing
% beta   = 1; alpha  = 5;  mu = 160; % good PSNR

% mu = max( alpha*norm(Wf(xbp),1), beta*linop_TV(Xbp) )/250;
solver = @(x0,z0,opts) solver_sBPDN_WW( A, alpha, W_tv, beta, W_wavelet, b, EPS, mu,x0, z0, opts);
solverF{testNo} = solver;

% k_end(testNo)   = 7;
% k_end(testNo)   = 1;
% maxits(testNo)  = 300;
%% for just tv
testNo  = 2;
string{testNo} = 'TV';
beta  = 1;
alpha = 0;
mu = 50;    % value used in paper
% mu = 30;  % optimal value of mu, after trying a few values
% mu = max( alpha*norm(Wf(xbp),1), beta*linop_TV(Xbp) )/250;

opts.normW2     = normTV^2;
solver = @(x0,z0,opts) solver_sBPDN_W( A, W_tv, b, EPS, mu, x0, z0, opts);
solverF{testNo} = solver;

% k_end(testNo)  = 7;
% k_end(testNo)  = 1;
% maxits(testNo) = 300;
%% for just analysis
testNo = 3;
string{testNo} = 'Analysis';
beta  = 0;
alpha = 1;
% mu = max( alpha*norm(Wf(xbp),1), beta*linop_TV(Xbp) )/500;
mu = 1; % value used in paper
opts.normW2     = normWavelet^2;
solver = @(x0,z0,opts) solver_sBPDN_W( A, W_wavelet, b, EPS, mu, x0, z0, opts);
solverF{testNo} = solver;
% k_end(testNo)   = 6;
% k_end(testNo)   = 1;
% maxits(testNo) = 300;
%% Solve with first-order method

for testNo = 1:3
    xPlug   = vec( Xbp );
    x       = xPlug;
    opts.maxits = maxits(testNo);
    opts.tol    = 5e-3;
    opts.stopCrit = 2;  % legacy option
    opts.errFcn = { @(f,dual,x) norm(vec(x)-vec(X_exact))/norm(X_exact,'fro') };
    lambda0     = [];
    tic
    % run continuation
    for k = 1:k_end(testNo)
        if k < k_end(testNo)
            opts.tol = opts.tol/2;
            opts.maxits = opts.maxits + 50;
        elseif k > 1
            opts.tol = 2e-4;
            opts.maxits = 500;
        end
        
        xPlug_old = x;
        solver   = solverF{testNo};
        [x,out,optsOut] = solver( vec(xPlug), lambda0, opts );
        lambda = cell(out.dual);
        xPlug_old = xPlug;
        xPlug = x + (k)/(k+3)*(x-xPlug_old);
        lambda0 = lambda;
    end
    X_recovered{testNo} = mat(x);
    toc


end

%%
figure(2);
splot = @(n) subplot(2,3,n);

splot(1);
imshow(X_exact);
title('Original');

splot(2);
imshow(Xbp);
title(sprintf('Noisy (PSNR %.1f dB)',...
    PSNR(Xbp)  ) );

splot(3);
imshow(X_oracle);
title(sprintf('Oracle wavelet threshold (PSNR %.1f dB)',...
    PSNR(X_oracle)  ) );

splot(4);
imshow(X_recovered{3});
title(sprintf('%s (PSNR %.1f dB)',...
    string{3}, PSNR(X_recovered{3})  ) );

splot(5);
imshow(X_recovered{2});
title(sprintf('%s (PSNR %.1f dB)',...
    string{2}, PSNR(X_recovered{2})  ) );

splot(6);
imshow(X_recovered{1});
title(sprintf('%s (PSNR %.1f dB)',...
    string{1}, PSNR(X_recovered{1})  ) );
% may not have best PSNR but visually looks very nice

