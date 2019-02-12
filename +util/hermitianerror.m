function nrm = hermitianerror(u0,u,p)
%HERMITIANERROR  Schatten p-norm error of  ||xx' - x0x0'||_p.
if nargin < 3
   p = 2;
end
[~,R] = qr([u0,u],0);
e = svd(R*blkdiag(-speye(size(u0,2)),speye(size(u,2)))*R');
nrm = norm(e, p);
end