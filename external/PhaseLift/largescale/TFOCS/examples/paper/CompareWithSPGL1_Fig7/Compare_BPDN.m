%{
    Compare the speed of SPGL1 to the Dual Conic formulation
using a Basis Pursuit and a Basis Pursuit De-Noising problem

To use SPGL1, need spgl1_modified.m and modified version of spgSetParms.m
These are included in the spgl1-1.7 directory, so please use this version
of SPGL1, not any version that may be pre-installed on your system.

Should take about 4 minutes to run this test

%}
%% Setup the path

addpath spgl1-1.7/
% Important: in the spgl1 directory, you may need to run spgsetup to compile the mex file

%% Test 5  (Haar wavelets of cameraman image) N = 65536 = 2^16
% Takes about 82 seconds to setup this test's reference solution
% and about 3.2 minutes to run all the tests
% This is figure 7 (b)
test=[];
test.name = 'test 5';
test.figName = 'Figure 7(b)';
test.seed = 9234;
test.SNR = 30;
test.sparse = false;
test.type = 'dct';
test.fileName = 'fig7b_reference';
test.xlim = [0,2000];
test.ylim = [1e-11,1];
cList   = [0,1];
TESTS{1} = test;
%% Test 6 Noiseless, N = 1024.
% Takes about 6 seconds to setup, and 18 seconds to run
% This is figure 7 (a)
test=[];
test.name = 'test 6';
test.figName = 'Figure 7(a)';
test.seed = 93825;
test.SNR = Inf;
test.DNR = 10;
test.sparse = true;
N = 2^14;
M = round(.5*N);
K = round(.25*M );
test.N = N; test.M = M; test.K = K;
test.type = 'DCT';
test.fileName = 'fig7a_reference';
test.xlim = [0,900];
test.ylim = [1e-14,1];
cList   = [0,2];
TESTS{2} = test;
%% Setup the tests:

for testCell = TESTS
  test = testCell{1};
  fprintf('Setting up exact solution for %s\n', test.name );
  if exist( [test.fileName,'.mat'], 'file')
      fprintf('Loading reference solution from %s\n', test.fileName );
      load(test.fileName);
      if ~isnumeric(Anew) && strcmpi(test.type,'dct')
          N = length(xNew); M = length(bNew);
          downsample = @(x) x(omega,:);
          SS.type = '()'; SS.subs{1} = omega; SS.subs{2} = ':';
          upsample = @(x) subsasgn( zeros(N,size(x,2)),SS,x);
          Af = @(x) downsample(idct(x.*D));
          At = @(x) D.*dct(upsample(x));
          Anew = {};
          Anew{1} = Af; Anew{2} = At;
      end
  else
      fprintf('Generating new reference solution\n');
      t2=tic;
      [xNew,Anew,bNew,normA2,EPS,xOrig,lambdaNew,D,omega] = setup_BPDN(test);
      toc(t2)
      save(test.fileName,'xNew','Anew','bNew','normA2','EPS',...
          'xOrig','lambdaNew','D','omega' );
  end
  Aff = Anew{1};
  Att = Anew{2};
  N = length(xNew); M = length(bNew);
  Af = @(x) multiplyA( Aff, Att, 'forward', x );
  At = @(x) multiplyA( Aff, Att, 'transpose', x );

  % === Check KKT conditions of reference solution
  % Lagrangian = ||x||_1 - (  <resid,lambda> + EPS*||lambda||_2 )
  resid = bNew - Aff(xNew);
  pdgap = resid'*lambdaNew + EPS*norm(lambdaNew);
  fprintf('Duality gap is\t\t%.2e\n', pdgap );
  fprintf('Primal infeasibility is\t%.2e\n', norm(resid) - EPS );
  % is 0 in subgradient of Lagrangian?
  T   = find(xNew);
  Tc  = find(~xNew);
  Atz = Att(lambdaNew);
  nrm1 = norm( sign(xNew(T)) + Atz(T), Inf ) ;
  nrm2 = norm( sign(xNew(T)) - Atz(T), Inf ) ; % some old conventions use -lamda
  fprintf('Stationarity: on T, subgradient condition violated by\t%.2e\n',...
      min( nrm1, nrm2 ) );
  fprintf('Stationarity: On T^c, subgradient condition violated by\t%.2e\n',...
      max( norm( Atz(Tc), Inf ) - 1, 0 ) );

  % ===  Run the tests and plot the results

  % CONTINUATION = false;   % no continuation
  % CONTINUATION = 1;       % "standard" continuation
  % CONTINUATION = 2;       % "accelerated" continuation (recommended)
  t1 = tic;
  for CONTINUATION = cList
      
      xPlug = Att(bNew)/normA2;
      alpha = norm(xPlug,1);   beta = norm(xPlug)^2/2;
      
      if strcmpi(test.name,'test 6')
          if CONTINUATION
              k_end = 4;
              mu0 = .1*(alpha/beta); % 0.2936
          else
              mu0 = .05*(alpha/beta); % 0.1468
          end
          tol = 1e-1*norm(bNew);
      end
      if strcmpi(test.name,'test 5')
          if CONTINUATION
              k_end = 30;
              mu0 = .1*(alpha/beta);
          else
              mu0 = .001*(alpha/beta);
          end
          tol = 1e-1*norm(bNew);
      end

      
      er        = @(x) norm(x-xNew)/norm(xNew);
      er_plug   = @(x) norm(x-xNew)/norm(xNew);
      
      % In these error computations, do NOT call multiplyA since we do NOT want
      % to count these.
      oExact    = norm(xNew,1);
      errp  = @(f,l,xk) norm( xk - xNew, 2 ) / norm( xNew,2 );
      
      globalErrors = @(f,l,xk) ComputeErrors( xNew, oExact, xk, f );
      multiplyA();  % zero it out.
      
      opts = [];
      opts.errFcn     = {errp, globalErrors};
      opts.Lexact     = normA2/mu0;
      x = zeros(N,1);
      xOld = x;
      x0 = x;
      opts.stopCrit = 2;  % legacy option
      lambda0  = [];
      
      
      opts.tol = tol;
      if CONTINUATION
          opts.printEvery = 100;
      else
          k_end = 1;  % no continuation
          opts.printEvery = 50;
          opts.maxits     = 2000;
      end
      
      continuation_locations = [];
      for k = 1:k_end
          if CONTINUATION == 2
              x0 = x + (k/(k+3))*(x-xOld); % accelerated continuation
              xOld = x;
          else
              x0 = x;
          end
          
          if k == k_end % at the final outer iteration, ask for more precision:
              opts.tol = 1e-11;
          else
              opts.tol = max( opts.tol/2, tol*1e-8 );
          end
          
          if EPS
              [x,out]=solver_sBPDN( linop_handles([M,N],Af,At), bNew, EPS,mu0,x0,lambda0,opts );
          else
              [x,out]=solver_sBP( linop_handles([M,N],Af,At), bNew,mu0,x0,lambda0,opts );
          end
          %     [x,lambda,out]=solve_BPDN( {Af,At}, bNew, mu0, EPS, opts );
          lambda = out.dual;
          lambda0 = lambda;
          
          eh = multiplyA('noreset');  % record error history
          continuation_locations = [continuation_locations; eh(end,1) ];
      end
      
      error_history = multiplyA();  % record error history
      
      ER{1 + CONTINUATION } = error_history;
      CONT_LOC{1 + CONTINUATION } = continuation_locations;
      
  end
  
  % --- solve via SPGL1
  Af = @(x) multiplyA( Aff, Att, 'forward', x );
  At = @(x) multiplyA( Aff, Att, 'transpose', x );
  options = [];
  options.errFcn = @(xk,f) ComputeErrors( xNew, oExact, xk, f );
  options.verbosity = 1;
  if strcmpi(test.name,'test 6')
      options.optTol = 1e-10;
  elseif strcmpi(test.name,'test 5')
      options.optTol = 1e-16;
  end
  multiplyA(); % zero it out.
  AA = @(x,mode) OperatorSGPL1(x,mode,Af,At);
  [x,r,g,info] = spgl1_modified( AA, bNew, 0, EPS, [], options );
  error_history = multiplyA();
  ER{4} = error_history;

  % ---- plot error_history
  % May not be smooth, due to backtracking, etc.
  figure; 
  handles = [];
  % list = setdiff( 1:length(ER), [3] );
  list = [1+cList,4];
  for i= list
      error_history = ER{i};
      iter = 1:length(error_history);
      h = semilogy( iter, error_history(:,2 ),'linewidth',2 );
      if i == 2 || i == 3
          % for continuation, add some extra plots
          clr = get(h,'color');
          continuation_locations = CONT_LOC{i};
          semilogy( iter(continuation_locations),...
              error_history(continuation_locations,2),...
              'og','color',clr,'markersize',10);
      end
      if i == 4  % SPGL1
          set(h,'linestyle','--');
      end
      handles = [handles,h];
      hold all
  end
  xlim(test.xlim);
  ylim(test.ylim);
  
  lString = {'AT w/o continuation','AT w/ continuation',...
      'AT w/ accel. continuation','SPGL1'};
  lh=legend(handles,lString(list));
  set(lh,'fontsize',20 );
  
  xlabel('number of calls to $A$ and $A^T$','interpreter','latex',...
      'fontsize',20);
  ylabel('relative error','fontsize',20);
  set(gca,'fontsize',16);
  title(test.figName);
  
  toc(t1);
  
end  % end of loop over tests

%{
    Error column legend:
column 1: iteration count
column 2: norm(x-xTrue)/norm(xTrue);
column 3: norm(x-xTrue,'inf');
column 4: length(T1)+length(T2)-2*length(intersect(T1,T2) );
column 5: abs( norm(x,1) -fTrue)/abs(fTrue);
%}
