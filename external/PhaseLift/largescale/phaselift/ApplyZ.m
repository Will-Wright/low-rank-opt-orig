function z = ApplyZ(v,z_old,Af,Pattern,g_Ay,stepsize,lambda,W)
%  ApplyZ computes Zv fast, used for eigs
%  here Z = X = stepsize*(A^* g_Ay + lambda *W*v
%  z_old, g_Ay, W need to be in factored form
%

% Author:     Thomas Strohmer 
%


Xv = z_old*(z_old'*v);


AtgAXv = AdjGradient(g_Ay,v,Pattern);

if iscell(W) == 0; % case W = I
  z = Xv - stepsize*(AtgAXv + lambda*v);
else % case W = (X + epsilon I)^{-1}
  epsilon = W{4};
  z = Xv - stepsize*(AtgAXv + lambda/epsilon*(v - W{1}*(W{2}'*v)));
end





% end of file
