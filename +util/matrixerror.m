function nrm = matrixerror(u0,v0,u,v,p)
%MATRIXERROR  Schatten p-norm error of  ||u0*v0' - u*v'||_p.
if nargin < 5
   p = 2;
end
[~,Ru] = qr([u0,u],0); [~,Rv] = qr([v0,v],0);
e = svd(Ru*blkdiag(-speye(size(u0,2)),speye(size(u,2)))*Rv');
nrm = norm(e,p);
end
