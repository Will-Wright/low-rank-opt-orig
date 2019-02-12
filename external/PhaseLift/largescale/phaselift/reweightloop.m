%% Reweighting (if needed) 
clear W
error_flag = 1;
lambda = 0.05;
if exist('maxloop','var') == 0;
  maxloop = 20;
end
iter = 0;
w_tol = 1e-3;
eigtol = 1e-7;
normxxt = sqrt(norm(y(:,1),1));  % this is the same as norm(x*x','fro')
                           % assuming the modulation operator is unitary

while error_flag
  iter = iter + 1;
  epsilon = normxxt/10;  % epsilon = norm(x*x','fro')/10
%  epsilon = 2;
  
  sqeig = normmat(X,0);
  tt = sqeig > eigtol;
  sqeig = sqeig(tt);
  % reduce epsilon in case it is too large, which can happen if
  % the first modulation is not unitary
  RNK = length(sqeig);
  if length(sqeig) == 1;    
     disp('Rank-one solution found - further reweighting not possible!');
     break; 
  end

  U = bsxfun(@rdivide,X(:,tt),sqeig);
  epsilon = min(epsilon,max(sqeig)/10);
  sqeig = 1./(sqeig + epsilon);
  W{1} = bsxfun(@times,U,sqrt(sqeig));
  W{2} = U;
  W{3} = sqeig;
  W{4} = epsilon;

  pars{1} = Pattern;
  pars{2} = lambda;
  pars{3} = epsilon;
  pars{4} = RNK;

  obj = @(U) AugNegLogLikelihood(U,Af,y,type,lambda,W);

  fprintf('\n\nReweighting\n');
  fprintf('');

  tic, [X,out,opts] = tfocs(obj, proj, pars, X, opts); toc

  % How well we do?
  err(iter) =  phaserr(x0,X(:,1));
  fprintf('\n\nHow well do we do?\n');
  fprintf('Approximate relative error = %9.4e\n', err(iter));
  svd(X)'
  res(iter) = norm(Af(X(:,1))-y)/norm(y)


  if err(iter) < w_tol; error_flag = 0; disp('stop'); end
  if iter > maxloop; break; end
end
