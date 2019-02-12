%% Setup problem size 
n = 50;              % signal is n x 1 (possibly complex valued) 
n_illumin = 6;       % number of structured illumination patterns

% Set signal and illumination patterns
x0      = randn(n,1) + 1i*randn(n,1);
Pattern = randn(n,n_illumin) + 1i*randn(n,n_illumin);

% Generate noiseless observations 
Af = @(X) StructuredIllumination(X,Pattern);
At = @(X) AdjStructuredIllumination(X,Pattern);
y = Af(x0*x0');
sz = { [n,n], [n,n_illumin] };

% Set up objective function and projection routine    
lambda = 0.1; W = eye(n);
obj = { @smooth_quad, smooth_linear( lambda, 0 ) };
linop = { linop_handles(sz,Af,At,'c2c'), -y ; ...
          linop_dot(W),                   0 };
proj = @proj_psd;

% Set starting point 
Xinit = zeros(n);

% Set TFOCS options
opts.alg        = 'AT';
opts.maxIts     = 2000;
opts.tol        = 1e-6;
opts.printEvery = 50;
opts.maxmin     = 1; % Minimize

% Call TFOCS!
[X,out,opts] = tfocs(obj,linop,proj,Xinit, opts);
% How well we do
% Print objective value for the different solutions
fprintf('\n\nHow well do we do?\n');
fprintf('Relative error = %9.4e\n', norm(x0*x0'- X,'fro')/norm(x0*x0','fro'));

%% Reweighting (if needed) 
epsilon = 2;
W = inv(X + epsilon*eye(n));
linop = { ...
    linop_handles(sz,Af,At,'c2c'), -y ; ...
    linop_dot(W), [] }; 
opts.tol        = 1e-8;

[X,out,opts] = tfocs(obj,linop,proj, X, opts);
% How well we do
fprintf('\n\nHow well do we do?\n');
fprintf('Relative error = %9.4e\n', norm(x0*x0'-X,'fro')/norm(x0*x0','fro'));
