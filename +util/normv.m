function normvx = normv(x,p)
   % NORMV computes the p-norm of a vectorized input
   if nargin < 2
      p = 2;
   end
   normvx = norm(x(:),p);
end
