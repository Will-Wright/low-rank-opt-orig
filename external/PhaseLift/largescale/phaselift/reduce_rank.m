function y = reduce_rank( x, z, theta )

% y = reduce_rank( x, z, theta )
% Computes the best rank-n approximation to (1-theta)*x + theta*z,
% where n  = size(x,2) = size(z,2);

[m,n] = size(x);
RNK = n;

eigsopts.issym = true;
eigsopts.isreal = false;

[V,D] = eigs(@(v) ApplyFactor(v,x,z,theta),m,RNK,'LM',eigsopts);

D = max(diag(D),0);
tt = D > 0;
tt = 1:length(D);
y  = bsxfun(@times,V(:,tt),sqrt(D(tt,:))');



