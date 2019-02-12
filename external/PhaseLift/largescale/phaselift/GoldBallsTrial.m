function[err,OutputX] = GoldBallsTrial(n_illumin, lambda, epsilon_reweight,filename,SNR,Pat,ReweightIter,n,yprime,Pattern)

%%% loads a 256x256 image of gold balls. image is scaled according to SNR
%%% 'high' or 'low'. epsilon_reweight sets epsilon for the re-weighting
%%% steps, it is initially set to Inf. When Pat == 0, 0/1 mask is used,
%%% otherwise a randn mask is used. Results saved in filename.

              % signal is n x n (possibly complex valued)
     % number of structured illumination patterns

 %load image of gold balls
 load goldballs_512DownSample.mat
 x0 = y; clear y;
 x0 = x0(1:n,1:n);
 
 
 y=yprime;
% %x0      = 1000*(randn(n,n) + 1i*randn(n,n));
% 
 

%%% choose SNR value
if(strcmp(SNR,'low'))
    SNRval = 1/5;
    SNRval = 10;
    x0 = SNRval*x0;
elseif(strcmp(SNR,'high'))
    SNRval = 60; 
    SNRval = 1100;
    x0 = SNRval*x0;
else
    disp('Error. SNR level undefined')
end

lambda = lambda*norm(x0(:));
 epsilon_reweight= epsilon_reweight*norm(x0(:));

% Pattern = zeros(n,n,n_illumin);
% 
% if(Pat == 0)   %%% 0/1 masks
% 
% for j = 1:n_illumin
% block = ones(4,1); mask1 = sign(randn(n/4,1)); mask1 = kron(mask1,block);
% block = ones(4,1); mask2 = sign(randn(n/4,1)); mask2 = kron(mask2,block);
% Pattern(:,:,j)= kron(mask1,mask2');
% end
% Pattern(Pattern==-1) = 0;
% Pattern(:,:,1) = 1;
% 
% %Pattern = rand(n,n,n_illumin);
% %Pattern(Pattern>0.5) = 1;
% %Pattern(Pattern ~=1) = 0;
% %Pattern(:,:,1) = 1;
% 
% elseif(Pat==1) %%% randn masks
% Pattern = randn(n,n,n_illumin) + 1i*randn(n,n,n_illumin);
% 
% else           %%% linear combination of sines mask
%     Pattern1 = SinePattern(n,n_illumin);
%     Pattern2 = SinePattern(n,n_illumin);
%     Pattern = zeros(n,n,n_illumin);
%     for j = 1:n_illumin
%         Pattern(:,:,j) = kron(reshape(Pattern1(:,j),n,1),reshape(Pattern2(:,j),n,1)');
%     end  
% end

% Generate noiseless observations
Af = @(X) FactoredIllumination(X,Pattern);

y_noiseless = Af(x0(:));

%norm(y(:)-y_noiseless(:))/norm(y_noiseless(:))

%y = poissrnd(y_noiseless,n*n,n_illumin);

%norm(y(:)-y_noiseless(:))/norm(y_noiseless(:))



% %y_length = length(y_noiseless(:));
% %gamma = (( 10^(SNRval/10))*(y_length)/(sum(sqrt(y_noiseless(:)))))^2;
% %x0 = sqrt(gamma)*x0;
% %y_noiseless = gamma*y_noiseless;
% 
% %%% Add poisson noise
% y = poissrnd(y_noiseless,n*n,n_illumin);
% 
% -20*log10(norm(y(:)-y_noiseless(:))/norm(y_noiseless(:)))
% 
% if(strcmp(SNR,'low'))
%     save('yLowSNR','y');
% elseif(strcmp(SNR,'high'))
%     save('yHighSNR','y')
% end



%%% For testing, try the algorithm with noiseless data
%y = y_noiseless;

%% Select model

%type = 'Poisson';
 type = 'Gaussian';

%% Select rank of matrix for truncated psd projection
RNK = 10;

%%  Set up objective function and projection routine   
 W = 1; %lambda = 0.05;


% Set starting point
Xinit = ones(n^2,RNK)/sqrt(n^2*RNK)*norm(x0(:));

Xinit1 = randn(n^2,1); Xinit2 = Xinit1(1:RNK); Xinit = Xinit1*Xinit2'; 
Xinit = Xinit/norm(Xinit(:))*norm(x0(:));

% Set TFOCS options
opts.alg        = 'AT';
opts.maxIts     = 1000;
opts.tol        = 1e-6;
opts.printEvery = 50;
opts.maxmin     = 1; % Maximize

proj = [];
epsilon = Inf;  % setting "epsilon = Inf" means W = Identity

obj = @(U) AugNegLogLikelihood(U,Af,y,type,lambda,W);

% setting additional parameters for reconstruction algorithm
pars{1} = Pattern;
pars{2} = lambda;
pars{3} = epsilon;
pars{4} = RNK;

% Call TFOCS!
tic, [X,out,opts] = tfocs(obj, proj, pars, Xinit, opts); toc
OutputX = X(:,1);
xyz = phaserr(x0,X(:,1));
err = zeros(1,ReweightIter+1);
  err(1) = xyz(1);
% How well we do?
fprintf('\n\nHow well do we do?\n');
fprintf('Approximate relative error = %9.4e\n', phaserr(x0(:),X(:,1)));

%size(X)
%X(:,size(X,2)+1:size(X,1)) = 0;

%[V,D] = eig(X); X = V*D^(1/2);
%rec = X(:,1);

y_trace = 0;

norms_trace = zeros(size(X,2));

%for j = 1:size(X,2)
%y_trace = y_trace + Af(X(:,j));
%if((norm(Af(X(:,j))) > 0.01) && (j<30))
%    norm(Af(X(:,j)))/norm(y_noiseless)
%end

%norms_trace(j) = norm(Af(X(:,j)));

%end
    
999999999

%norm((y_trace-y_noiseless))/norm(y_noiseless)

%[V,D] = eig(X); Xfact = V*sqrt(D)*V';
%y_trace = Af(Xfact);

%measError = norm(y_trace(:)-y_noiseless)


%reweightloop;

%% Reweighting (if needed) 
clear W
large_error = 1;
maxloop = 5;
iter = 0;
w_tol = 1e-3;
eigtol = 1e-6;



for j = 1:ReweightIter
  iter = iter + 1;
  epsilon = epsilon_reweight; 
  %epsilon = 2;
  
  %lambda = lambda/2;
  sqeig = normmat(X,0);
  tt = sqeig > eigtol;
  sqeig = sqeig(tt);
  RNK = length(sqeig);
  if length(sqeig) == 1;  
     disp('Rank-one solution found - further reweighting not possible!');
    % break; 
  end

  U = bsxfun(@rdivide,X(:,tt),sqeig);
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
  Outputs(iter+1) = {X};
  
  
save('CheckThisShit','Outputs','x0','y_noiseless','y');

  % How well we do?
  xyz = phaserr(x0,X(:,1));
  err(iter+1) = xyz(1);
  OutputX = X(:,1);
  
  fprintf('\n\nHow well do we do?\n');
  fprintf('Approximate relative error = %9.4e\n', err(iter));
  

  %if err(iter) < w_tol; large_error = 0; disp('stop'); end
  %if iter > maxloop; disp('need to break'); break; end
end


%save(filename,'err')




