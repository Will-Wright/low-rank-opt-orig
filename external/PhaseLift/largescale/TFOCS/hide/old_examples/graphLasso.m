%{
Estimates a partially observed covariance matrix
by solvingn the "graph LASSO".

Let S be the empiricaly covariance matrix, then the MLE is

X = argmin <X,S> - log det(X)
    subject to X > 0

which is solvable in closed form (since S >= 0, and usually S > 0)
by X = inv(S)

We add a sparsity constraint on X and solve:

X = argmin <X,S> - log det(X) + lambda ||X||_1
    subject to X > 0
where ||X||_1 is component-wise.
We can't compute the prox to ||X||_1 and X > 0 simultaneously though...

Need to solve dual :-(


So, this is a bit more complicated....


Assume we have zero mean

%}
N = 50;     % dimensions
m = 5000;    % number of observations

% Not sure if this is a reasonable setup
SigmaInv    = sprand(N,N,.1/2);
SigmaInv    = SigmaInv*SigmaInv';   % symmetric, and positive semi-definite
SigmaInv    = SigmaInv - diag(diag(SigmaInv)) + 100*diag(ones(N,1));
fprintf('Sparsity of inv(Sigma) is %.2f%%\n', 100*(nnz(SigmaInv)/numel(SigmaInv)));
spy(SigmaInv);
Sigma       = inv(full(SigmaInv));
R           = chol(Sigma);
Z           = randn(m,N)*R;

S           = cov(Z);
invS        = inv(S);
fprintf('Error using S is %.2f\n', norm(S-Sigma)/norm(Sigma) );
fprintf('Error using inv(S) is %.2f\n', norm(invS-SigmaInv)/norm(SigmaInv) );

%% Solve in TFOCS
opts = [];
opts.printEvery     = 10;
opts.tol            = 1e-10;
opts.errFcn     = @(f,x) norm(x-SigmaInv)/norm(SigmaInv);
lambda  = 1;
X0      = inv(S);
% X = tfocs( {smooth_linear(S); smooth_logdet},{[];[]}, {prox_l1(lambda),
% proj_psd}, [], opts );
X = tfocs( {smooth_linear(S), smooth_logdet},{1;1},proj_psd, X0, opts );