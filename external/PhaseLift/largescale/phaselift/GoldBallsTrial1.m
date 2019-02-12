function[y,Pattern]  = GoldBallsTrial1(n_illumin, lambda, epsilon_reweight,filename,SNR,Pat,ReweightIter,n)

load goldballs_512DownSample.mat
x0 = y; clear y;
x0 = x0(1:n,1:n);
%x0      = 1000*(randn(n,n) + 1i*randn(n,n));



%%% choose SNR value
if(strcmp(SNR,'low'))
    %used to be 8.6
    SNRval = 1/5; 
    SNRval = 10;
    x0 = SNRval*x0;
elseif(strcmp(SNR,'high'))
    %used to be 11
    SNRval = 60; 
    SNRval = 1100;
    
    
    x0 = SNRval*x0;
else
    disp('Error. SNR level undefined')
end

lambda = lambda*norm(x0(:));
epsilon_reweight= epsilon_reweight*norm(x0(:));

Pattern = zeros(n,n,n_illumin);

if(Pat == 0)   %%% 0/1 masks

%for j = 1:n_illumin
%block = ones(4,1); mask1 = sign(randn(n/4,1)); mask1 = kron(mask1,block);
%block = ones(4,1); mask2 = sign(randn(n/4,1)); mask2 = kron(mask2,block);
%Pattern(:,:,j)= kron(mask1,mask2');
%end
%Pattern(Pattern==-1) = 0;
%Pattern(:,:,1) = 1;

Pattern = rand(n,n,n_illumin);
Pattern(Pattern>0.5) = 1;
Pattern(Pattern ~=1) = 0;
Pattern(:,:,1) = 1;

elseif(Pat==1) %%% randn masks
Pattern = randn(n,n,n_illumin) + 1i*randn(n,n,n_illumin);
Pattern(:,:,1) = 1;

else           %%% linear combination of sines mask
    Pattern1 = SinePattern(n,n_illumin);
    Pattern2 = SinePattern(n,n_illumin);
    Pattern = zeros(n,n,n_illumin);
    for j = 1:n_illumin
        Pattern(:,:,j) = kron(reshape(Pattern1(:,j),n,1),reshape(Pattern2(:,j),n,1)');
    end  
end

% Generate noiseless observations
Af = @(X) FactoredIllumination(X,Pattern);

y_noiseless = Af(x0(:));
%y_length = length(y_noiseless(:));
%gamma = (( 10^(SNRval/10))*(y_length)/(sum(sqrt(y_noiseless(:)))))^2;
%x0 = sqrt(gamma)*x0;
%y_noiseless = gamma*y_noiseless;

%%% Add poisson noise
%y = poissrnd(y_noiseless,n*n,n_illumin);
%y = normrnd(y_noiseless,n*n,n_illumin);

%%% Add gaussian noise

%y = y_noiseless;% + normrnd(0,norm(x0)/,size(y_noiseless));
%%% 32 masks y = y_noiseless + normrnd(0,norm(x0)/0.068,size(y_noiseless));
y = y_noiseless + normrnd(0,norm(x0)/0.048,size(y_noiseless));


%y = y_noiseless + normrnd(0,norm(x0)/0.6,size(y_noiseless));


display('SNR')
-20*log10(norm(y(:)-y_noiseless(:))/norm(y_noiseless(:)))

norm(y(:)-y_noiseless(:))/norm(y_noiseless(:))

