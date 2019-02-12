AddPathThomasTFOCS

dim = 6;
n_illumin = 8;
lambda = 1/500; %%% scale of lambda (to be multiplied by the norm of the signal)
epsilon_reweight = 1/10; %%% scale of epsilon (to be multiplied by the norm of the signal)
filename = '2B';
SNR = 'high';
Pat = 1; %%% Pat = 0 is for 0/1 masks, Pat = 1 is for randn masks, Pat = 2 is for sinusoid masks
ReweightIter = 5;
NumBisection = 7; 


[y,Pattern] = GoldBallsTrial1(n_illumin, lambda, epsilon_reweight,filename,SNR,Pat,ReweightIter,dim);

%f = @(l) GoldBallsTrial(n_illumin, l, epsilon_reweight,filename,SNR,Pat,ReweightIter,dim,y);


[BestMSE, BestX] = GoldBallsTrial(n_illumin, lambda, epsilon_reweight,filename,SNR,Pat,ReweightIter,dim,y,Pattern);

%[BestMSE,BestX] = FindBestLambda(f,ReweightIter,0,lambda,0,NumBisection,dim);

BestMSE

save('2BHighSNR8nLambda1TEST','BestMSE','BestX')

%%