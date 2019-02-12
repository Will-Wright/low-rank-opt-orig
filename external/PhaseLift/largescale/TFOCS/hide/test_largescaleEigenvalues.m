%{

    This tests largescale problems, where you may only
    want to calculate an approximate gradient.

    For now, this is very small scale, using CVX as a reference solution,
    but we make sure to use function handles and factored forms so that
    it is ready to scale

Stephen Becker, Nov 1 2010

%}

randn('state',324324);
rand( 'state',234324);

n = 16;
N = n^2;
M = round(N/2);

% our linear operator, which acts on X = WW^*
A = randn(M,N);
Af = @(W,Wt) A*