function y = innprodW(X,W)
% INNPRODW(W,X) - computes <X,W> where the matrices X and W are assumed to
%                 be given in factorized form
%
%
% Input:  
%             
% Output:    
%

% Author:     T.Strohmer (strohmer@math.ucdavis.edu)
%

if iscell(W) == 1;
  epsilon = W{4};
  y = 1/epsilon*(norm(X,'fro')^2 - norm(X'*W{1},'fro')^2);
else % case W = (X + epsilon*I)^{-1}
  y = norm(X,'fro')^2;
end





% end of file
