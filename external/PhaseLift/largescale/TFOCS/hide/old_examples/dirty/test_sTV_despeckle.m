%{
    Despeckle an image:

    min_{x,y} alpha*||x||_TV + beta*||y||_1
s.t.
    || A(x+y) - b || <= eps

The solvers solve a regularized version

see also test_sTV_largescale.m and test_sTV_Analysis_largescale.m

This version has no reference solution.

NOT WORKING, DO NOT INCLUDE!!!!!!

%}

% Before running this, please add the TFOCS base directory to your path


% Generate a new problem

randn('state',245);
rand('state',245);

% image = 'phantom';
image = 'bugsbunny';

switch image
    case 'phantom'
        n = 256;
        n1 = n;
        n2 = n;
        N = n1*n2;
        n = max(n1,n2);
        x = phantom(n);
        x = x(1:n1,1:n2);
        
        p = .01; % percent of corrupted pixels
    case 'bugsbunny'
        % Image is from wikimedia commons:
        % http://commons.wikimedia.org/wiki/File:Falling_hare_bugs.jpg
        load BugsBunny
        x = bugs;
        [n1,n2] = size(x);
        N = n1*n2;
        
        p = .01; % percent of corrupted pixels
end
x_original = x;

mat = @(x) reshape(x,n1,n2);
vec = @(x) x(:);

M = N;
% Our measurement is "X+Y", but we treat this as a single variable [X;Y]
% so "A" must split it apart:

A_matrix    = [speye(N), speye(N)];
A           = linop_matrix(A_matrix);
normA       = sqrt(2);
   
b_original = vec(x_original);

% Add noise via a random mask:
x_temp  = awgn( x_original, 0, 'measured'); % adds a lot of noise
index   = randperm(N);
index   = index( 1 : round(p*N) );
x_noisy = x_original;
x_noisy(index) = x_temp(index);

b   = vec(x_noisy);
EPS = .9*norm(b-b_original);
figure(2); imshow( [x_original, x_noisy] );

%% Setup "W"
% We break the primal variable into two parts: [x;y]
% So "W" just picks out y.
W_y = linop_matrix( [sparse(N,N), speye(N) ] ); % pick out Y variable
normW_y     = linop_test(W_wavelet);

%% Setup TV
W_x = linop_matrix( [speye(N), sparse(N,N) ] ); % pick out X variable
W_tv= linop_TV( [n1,n2] );
W_tv= linop_compose( W_tv, W_x );
normW_x     = linop_test(W_tv);

%% Call the TFOCS solver

x_ref   = x_original;

norm_x_ref      = norm(x_ref,'fro');
er_ref          = @(x) norm(vec(x)-vec(x_ref))/norm_x_ref;
er              = er_ref;  
opts = [];
opts.restart    = 1000;
% opts.errFcn     = @(f,dual,primal) er(primal(1:N)+primal(N+1:end));
opts.errFcn     = @(f,dual,primal) er(primal(1:N) );
opts.maxIts     = 500;
opts.printEvery = 20;


% x0 = [ vec(x_noisy); zeros(N,1) ];
x0 = [vec(x_noisy)*.9; .1*vec(x_noisy) ];
z0  = [];   % we don't have a good guess for the dual

% A   = W_x;
% normA   = normW_x;

opts.normW12     = normW_x^2;
opts.normW22     = normW_y^2;
opts.normA2      = normA^2;


% Make sure we didn't do anything bad
% linop_test( A );
% linop_test( W );

contOpts            = [];
contOpts.maxIts     = 2;
contOpts.betaTol    = 2;

mu      = 1;
alpha   = 1e-10;
beta    = 1e-10;

opts.tol        = 1e-8*min([alpha,beta]);

disp('============ DESPECKLING =============');
epsilon = 0;
% epsilon = EPS;

if epsilon
    proxA   = prox_l2(epsilon);
else
    proxA   = proj_Rn;
end

nA = linop_compose(A,-1);

W1 = W_x;
W2 = W_y;
nW2 = linop_compose(W2,-1);
proxScale1  = sqrt( opts.normW12 / opts.normA2 );
proxScale2  = sqrt( opts.normW22 / opts.normA2 );
prox        = { proxA, ...
                proj_linf(proxScale2*beta) };%,...
%                 proj_linf(proxScale2*beta) };%, ...
% prox        = { proxA, ...
%                 proj_linf(proxScale1*alpha),...
%                 proj_linf(proxScale2*beta) };%, ...
%                 proj_Rn, ...
%                 proj_Rn };
W1          = linop_compose( W1, 1 / proxScale1 );
W2          = linop_compose( W2, 1 / proxScale2 );
% solver = @(mu,x0,z0,opts)tfocs_SCD( [], { A, -b; W1, 0; W2, 0}, ... %;
% ...
solver = @(mu,x0,z0,opts)tfocs_SCD( [], { A, -b; W2, 0}, ... %; ...
    prox, mu, x0, z0, opts );
%     W2, 0*ones(N,1); nW2, ones(N,1)  }, ...
    

% solver = @(mu,x0,z0,opts) solver_sBPDN_WW( A, alpha, W_x, beta, W_y, b, epsilon, mu,x0, z0, opts);

[ x, out, optsOut ] = continuation(solver,mu,vec(x0),z0, opts,contOpts);
xx = x(1:N);
xy = x(N+1:end);
% z   = xx + xy;
z = xx;
X_despeckled = mat( z );
X_speckled   = mat( xy );

[norm(xx),norm(xy),nnz(thresh(xy))/numel(xy) ]

%
figure(1); clf;
splot = @(n) subplot(2,2,n);

splot(1);
imshow(x_original);
title(sprintf('Noiseless image, PSNR %.1f dB', Inf ));

splot(2);
imshow(x_noisy);
title(sprintf('Noisy image, PSNR %.1f dB', PSNR(x_noisy) ));

% splot(3);
% imshow(X_hat);
% title(sprintf('Oracle wavelet thresholding, PSNR %.1f dB', PSNR(X_hat) ));

% splot(4);
% imshow(X_wavelets);
% title(sprintf('Wavelet regularization, PSNR %.1f dB', PSNR(X_wavelets) ));
% 
% splot(5);
% imshow(X_tv);
% title(sprintf('TV regularization, PSNR %.1f dB', PSNR(X_tv) ));

splot(3);
imshow(X_despeckled);
title(sprintf('Despeckling, PSNR %.1f dB', PSNR(X_despeckled) ));

splot(4);
imshow(X_speckled);
title(sprintf('Speckles'));

%%
% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
