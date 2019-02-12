%{
    Example of the Random Modulator Pre-Integrator (RMPI)
    with a radar pulse and continuous wave interferer

    This was used to generate figure 11 in the paper
    Stephen Becker

    I used the communication toolbox, but it shouldn't be necessary
        (below I have provided work-around functions).
    Also used the signal processing toolbox.  I think this is necessary,
        since recreating the hamming window and periodogram functions
        takes more than trivial work.

    Also relies on some helper files, which are included in this directory
    these files are:
        PsiWff and PsiTransposeWff  (the Gabor dictionary code)
        pulseTrain      (for setting up a sample pulse window )


    Should take about 4 or 5 minutes to run BP and a similar time
    to run Dantzig, so about 10 minutes to run this whole file.

%}

%-- test for toolboxes
% clear awgn randsrc
%if ~license('test','communication_toolbox')
    awgn = @(x,snr,ignore) x + ...
        10^( (10*log10(sum(abs(x(:)).^2)/length(x(:))) - snr)/20 )*randn(size(x));
    randsrc = @(varargin) sign( rand(varargin{:}) -.5 );
%end

%-- setup signal
factor = 4;
N   = factor*2048;
FS  = 5e9;  % Nyquist rate
BW  = FS/2; % bandwidth if we sampled at Nyquist rate
dT  = 1/FS;
KHz = 1e3; MHz = 1e6; GHz = 1e9; ms = 1e-6; ns = 1e-9; ps = 1e-12;
T   = N*dT;
t   = (0:N-1)*dT;


randn('state',324324); rand('state',324234);
% -- Add a large pu;se "x1"
f1  = rand*BW;  % frequency
a1  = 1;        % amplitude
p1  = rand*2*pi;% phase
% -- window by a trapezoidal window
tArrive = factor*40*ns;
tRise   = factor*20*ns;  tFall = tRise; 
tDur    = factor*300*ns;
tRep    = 0;
% -- get a trapezoidal window:
envelope1 = pulseTrain( t, tArrive, -tRise, tDur, tFall, tRep );
% -- this is unrealistic, so smooth it a bit:
f_cutoff = .004*FS;  % location of the first zero in frequency
n = find( t < 2/2/f_cutoff ); n = n(end);
kernel = ones(n,1); kernel = kernel / sum(kernel);
envelope1 = conv( envelope1, kernel, 'same' )';
x1  = a1*envelope1.*(sin(2*pi*f1*t + p1 ))';  % signal

% Add a small pulse "x2"
f2  = rand*BW;
a2  = 1e-3*a1;
p2  = rand*2*pi;
s   = a2*sin(2*pi*f2*t + p2);
% -- window by a trapezoidal window
tArrive = factor*45*ns;
tRise   = factor*20*ns;  tFall = tRise; 
tDur    = factor*300*ns;
tRep    = 0;
% -- get a trapezoidal window:
envelope2 = pulseTrain( t, tArrive, -tRise, tDur, tFall, tRep );
% -- this is unrealistic, so smooth it a bit:
f_cutoff = .004*FS;  % location of the first zero in frequency
n = find( t < 2/2/f_cutoff ); n = n(end);
kernel = ones(n,1); kernel = kernel / sum(kernel);
envelope2 = conv( envelope2, kernel, 'same' )';
x2 = envelope2.*s';

% -- combine big and small signals:
x_exact = (x1 + x2);

% -- plot the pulse windows
figure(1); clf;
xDemod = a1*envelope1;
plot( t/ns, xDemod,'linewidth',2 );
hold all
xDemod = 100*a2*envelope2;
plot( t/ns, xDemod,'linewidth',2 );
xlabel('time in ns','fontsize',20);
ylabel('amplitude of envelope','fontsize',20);

h = legend('large pulse','small pulse, scale exaggerated by 100');
set(h,'fontsize',20);
ylim([-.1,1.1])



%% measurements
randn('state',320424); rand('state',3247825);
F_ADC       = 50*MHz;
nChannels   = 8;    % so, equivalent to a 400 MHz sampler
nSamples    = round(FS/F_ADC);
nPeriods    = floor( N/nSamples);
M           = nPeriods*nChannels;
Phi         = zeros(M,N);
for j = 0:nPeriods-1
    Phi(   j*nChannels + (1:nChannels) , j*nSamples + (1:nSamples)  ) = ...
        randsrc(nChannels,nSamples);
end
% add a bit to the final period
Phi( j*nChannels + (1:nChannels), nPeriods*nSamples:N ) = ...
    randsrc( nChannels, N - nPeriods*nSamples+1 );
normA2  = norm(Phi*Phi');
Phi = Phi/sqrt(normA2);

A   = @(x) Phi*x;
At  = @(y) Phi'*y;
b_exact = A(x_exact);
SNR = 60;   % in dB
x_noisy = awgn(x_exact, SNR, 'measured' );
z = x_exact - x_noisy;
b       = A(x_noisy);
EPS = norm(b-b_exact);
normA2  = norm(Phi*Phi');
fprintf('EPS is %.3f\n', EPS )

win = hamming(N);
clf; periodogram( x_noisy,win, [], FS );

rms = @(x) norm(x)/sqrt(length(x));
snr = @(sig,noise) 20*log10( rms(sig)/rms(noise) );

fprintf('SNR of signal is %.1f dB; SNR of big tone is %.1f dB; SNR of small tone is %.1f dB\n',...
    snr(x_exact,z), snr(x1,z), snr(x2,z) );

sigma = std( b - b_exact );


%% setup dictionary
% The PsiWFF and PsiTransposeWFF code is a Gabor frame
% (i.e. a short-time Fourier transform)
% written by Peter Stobbe
% 
% PsiWFF is the synthesis operator, and acts on coefficients
% PsiTransposeWFF is the analysis operator, and acts on signals
gMax = 0;
gLevels = 5;
tRedundancy = 1;
fRedundancy = 1;
gWindow = 'isine';

lString = {};
figure(2); clf;

% gabor.gMax = gMax; gabor.gLevels = gLevels;
% gabor.tRedundancy = tRedundancy; gabor.fRedundancy = fRedundancy;
% param.gabor = gabor;

logN = log2(N);
psi_A = @(y) PsiWFF(y,gWindow,logN,logN-gLevels,logN+gMax,tRedundancy,fRedundancy);
psiT_A = @(x) PsiTransposeWFF(x,gWindow,logN,logN-gLevels,logN+gMax,tRedundancy,fRedundancy);

N_Gabor = length( psiT_A( ones(N,1) ) );
% param.N_Gabor = N_Gabor;

% DO we need this counter?
% psi = @(y) counter( psi_A, y);
% psiT= @(x) counter( psiT_A,x);
psi  = psi_A;
psiT = psiT_A;

x  = x_exact;
y  = psiT(x);
x2 = psi(y);
d  = length(y);
% The frame is tight, so the psuedo-inverse is just the transpose:
fprintf('Error in PSI(PSI^* x ) - x is %e; N = %d, d = %d, d/N = %.1f\n', ...
    norm(x-x2),N,d,d/N);

normW2 = 1;

% look at coefficients
fcn = @(x) sort(abs(x),'descend');
c = fcn( psiT(x) );
figure(2); 
%clf;
semilogy(c)
hold all

% estimate a power-law decay rate, x_{k} = a*k^p
% offset = round(.2*d);
offset = 1;
cBulk = c( offset:round(.9*d) );
K = length(cBulk);
k = (1:K);
a = [ones(K,1), log(1:K)'] \ log(cBulk);
p = a(2);
% hold all
% semilogy( offset-1+ k, exp(a(1))*k.^p );
lString{end+1} = sprintf('p = %.3f',p );
fprintf('  and p is %.3f\n', p );

% end
legend(lString)
%% for Dantzig, need to find D and delta
% DANTZIG = true;

for DANTZIG = 0:1

  if DANTZIG
    alpha   = .05;
    AtA     = Phi'*Phi;
    Anorms  = sqrt( sum(AtA.^2) )';
    nTrials = min(4*N,400);
    w       = randn(N,nTrials);
    Aw      = Phi'*( Phi*w);
    supAtz  = sort(max(Aw ./Anorms(:,ones(1,nTrials))));
    thresh  = supAtz(round(nTrials*(1-alpha)));  % empirical
    sigma   = std( x_exact - x_noisy);
    delta   = thresh*sigma*Anorms;  % a vector, for Dantzig solvers
    normA2  = norm(Phi*Phi');
    if all(delta > 1e-10)
        delta0  = mean(delta)/normA2;
        D       = delta0./delta; % watch out for division by 0
    else
        D       = 1;
        delta0  = 0;
    end
    normA2 = norm(D,Inf)^2;
    clear alpha w thresh nTrials
  end

  %--  RECONSTRUCT via dual conic smoothing
  W  = psiT;
  Wt = psi;
  opts        = [];
  opts.maxits = 2500;
  opts.maxmin = -1;
  opts.saddle = 1;
  opts.printEvery = 150;
  opts.tol    = 1e-1;
  opts.errFcn = { @(l,f,x) norm(x-x_exact)/norm(x_exact) };
  if DANTZIG
      delta = 0.5*delta0;
      opts.errFcn{end+1} = @(l,f,x) norm(D.*(At(A(x)-b)),Inf)-delta;
  else
      opts.errFcn{end+1} = @(l,f,x) norm(A(x)-b)-EPS;
  end
  
  solver  = @solver_AT;
  xx      = Phi\b;
  alpha   = norm( W(xx), 1 );
  beta    = norm( xx )^2/2;
  mu      = .1*(alpha/beta);
  
  scaleW  = sqrt(normA2/normW2);
  MU      = scaleW*mu;
  opts.L0 = 2/MU*normA2;
  opts.stopCrit   = 2; % legacy
  
  lambda0     = [];
  
  xPlug   = zeros(N,1);
  xOld    = xPlug;
  xOld_rw = xOld;
  
  % reweighting loop
  n_rw = 8;
  time=[]; errs = []; iters = [];
  
  % Here are the relevant parameters
  tol_rw              = 1e-2;
  tol_continuation    = 1e-3;
  cutoff_rw           = 95;  % of 100
  tol_increase        = 1.5;


  tt = clock;
  for rw = 1:n_rw  % reweighting loop
    
    WW = linop_handles([N_Gabor,N],W,Wt );
    fprintf('\n\n-- Reweighting %d of %d --\n\n', rw, n_rw );

    
    opts.tol    = 1e-1;
    k_end   = 12;
    for k = 1:k_end  % continuation loop       
        
        tic;
        opts.normA2     = normA2;
        opts.normW2     = normW2;
        if DANTZIG
            [x,out,optsOut] = solver_sDantzig_W({Phi,D},WW,b,delta,mu,xPlug,lambda0,opts );
        else
            [x,out,optsOut] = solver_sBPDN_W(Phi,WW,b,EPS,mu,xPlug,lambda0,opts );
        end
        time(rw,k) = toc;
        errs(rw,k) = out.err(end,1);
        iters(rw,k) = out.niter;
        
        lambda = out.dual;
        fprintf('Relative error is %.2e\n', norm(x-x_exact)/norm(x_exact) );
        
        % continuation update
        xPlug   = x + k/(k+3)*( x - xOld );
        lambda0 = lambda;
        
        opts.tol = opts.tol/tol_increase;
        if k > 1
          change = norm( x - xOld )/norm(x);
          fprintf('After continuation %d, rel. change in x is %.2e\n',k,change);
          if change < tol_continuation, disp('Change is small, so exiting from continuation'); break; end
        end
        xOld    = x;
        
        
    end
    
    % update weights
    coeff   = psiT(x);
    cSort   = sort(abs(coeff),'descend');
    ind     = [find( cumsum(cSort.^2)/sum(cSort.^2) > (cutoff_rw  )/100 );length(cSort)];
    fprintf('cutoff index is %d of %d total, i.e. %.4f%%\n', ind(1),length(cSort),...
        100*ind(1)/length(cSort) );
    cutoff  = cSort(ind(1));
    weights = 1./( abs(coeff) + cutoff );
    W       = @(x) weights.*psiT(x);
    Wt      = @(y) psi(weights.*y);
    
    change  = norm( x - xOld_rw)/norm(x);
    fprintf('After reweighting %d, rel. change in x is %.2e\n',rw, change );
    if change < tol_rw, disp('Change is small, so exiting from reweighting'); break; end
    xOld_rw = x;
    xPlug   = x;
    
    
  end
  if DANTZIG, type = 'dantzig'; else type = 'bpdn'; end
  %save(sprintf('rmpi_results_%s',type),'time','errs','iters');
  fprintf('TOTAL ELAPSED TIME IS %.1f MINUTES\n', etime(clock,tt)/60 );

  
  %-- Make plot
  figure(1+DANTZIG); clf;
  win = hamming(N);
  NN = 1*N;
  
  [Pxx,F] = periodogram( x_exact, win, NN, FS );
  Pxx = db(Pxx,'power');
  F = F/GHz;
  plot( F, Pxx,'k','linewidth',2);
  hold all
  
  [Pxx,F] = periodogram( x, win, NN, FS );
  Pxx = db(Pxx,'power');
  F = F/GHz;
  plot( F, Pxx,'r.','markersize',20 );
  
  xlabel('Frequency (GHz)','fontsize',20);
  ylabel('Power/frequency (dB/Hz)','fontsize',20);
  
  ylim([-220,-50])
  
  h = legend('Original signal','Recovered signal');
  set(h,'fontsize',20 );
  
  % figure(3); clf;
  % plot( x )

end
%% Display data for tables
for DANTZIG = 0:1
  if DANTZIG, type = 'dantzig'; else type = 'bpdn'; end
  load(sprintf('rmpi_results_%s',type))
  %--  make a table showing time and # of iterations
  if DANTZIG, type = 'dantzig'; else type = 'bpdn'; end
  rw = size(time,1);
  
  fprintf('For %s\n', type);
  for rw_i = 1:rw
      nCont = length( find(iters(rw_i,:) ) );
      nn = sum(iters(rw_i,: ));
      tm = sum(time(rw_i,:));
      e  = errs(rw_i,nCont);
      fprintf('Reweight %d, %d continuation steps, %d total iters, %.1f seconds, %.2e final error\n',...
          rw_i, nCont, nn, tm, e );
  end
end


