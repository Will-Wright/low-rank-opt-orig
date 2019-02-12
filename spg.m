function [x, stats] = spg(FG, P, x, opts)
% SPG  Spectral projected gradient.

   % Parse options.
   if nargin == 3 || ~isstruct(opts)
      opts = [];
   end
   if isempty(P)
      P = @(X)X;
   end
   if ~isfield(opts,'pgtol') || isempty(opts.pgtol)
      opts.pgtol = 1e-6;
   end
   if ~isfield(opts,'maxit') || isempty(opts.maxit)
      opts.maxit = 2^14;
   end
   if ~isfield(opts,'useclassicalgradients') || isempty(opts.useclassicalgradients)
      opts.useclassicalgradients = true;
   end
   if ~isfield(opts,'verbose') || isempty(opts.verbose)
      opts.verbose = false;
   end
   if ~isfield(opts,'ftol') || isempty(opts.ftol)
      opts.ftol = -inf;
   end
   
   normv = @(v)norm(v(:));
   stats.time = tic;
   y = P(x); [fy,gy] = FG(y); x = reshape(x,size(gy)); y = reshape(y,size(gy)); ny = x-y; ay = 1;
   stats.xInit =  y; stats.xOut = stats.xInit;
   stats.fInit = fy; stats.fOut = stats.fInit;
   stats.gInit = gy; stats.gOut = stats.gInit;
   stats.nInit = ny; stats.nOut = stats.nInit;
   stats.aInit = ay; stats.aOut = stats.aInit;
                     stats.kOut = 0;
   stats.stepsEstimated = [];
   stats.stepsFound     = [];
   stats.pgNorms        = normv(y-P(y-gy));
%  pgTol = opts.pgtol * max(1,stats.pgNorms(end));
   pgTol = opts.pgtol * (1+stats.pgNorms(end));
   fTol = opts.ftol;
   stats.countFG = 1;
   stats.countP  = 2;
   if ~isfield(opts,'ls') || ~isstruct(opts.ls)
      opts.ls = [];
   end
   stats.functions     = fy;
   stats.ls.functions  = [];
   stats.ls.iterations = [];

   if opts.verbose
      fprintf('%4s %4s %10s %10s %10s\n','iter','lsit','fobj','normpg','feas');
   end

   k = 0;
   while k < opts.maxit
      [x,fx,gx,nx,ax] = deal(y,fy,gy,ny,ay);
      stats.stepsEstimated = [stats.stepsEstimated;ax];
      [y,fy,gy,ny,ax,opts.ls] = linesearch(FG,P,x,fx,gx,nx,ax,opts.ls);
      stats.pgNorms       = [stats.pgNorms;normv(y-P(y-gy))];
      stats.functions     = [stats.functions;fy];
      stats.ls.functions  = [stats.ls.functions ;opts.ls.functions ];
      stats.ls.iterations = [stats.ls.iterations;opts.ls.iterations];
      stats.stepsFound = [stats.stepsFound;ax];
      stats.countFG = stats.countFG+opts.ls.iterations;
      stats.countP  = stats.countP +opts.ls.iterations;
      if opts.useclassicalgradients
         ay = estimatestep(y-x,gy-gx,opts.ls);
      else
         ay = estimatestep(y-x,gy+ny-gx-nx,opts.ls);
      end
      if isnan(fy) || ~isfinite(fy)
         break;
      end
      if fy < stats.fOut
         stats.xOut =  y;
         stats.fOut = fy;
         stats.gOut = gy;
         stats.nOut = ny;
         stats.aOut = ay;
         stats.kOut = k;
      end
      
      % Logging.
      if opts.verbose
         fprintf('%4d %4d %10.2e %10.2e %10.2e\n',k,opts.ls.iterations,fy,stats.pgNorms(end)/stats.pgNorms(1),normv(y-P(y)));
      end
      
      k = k + 1;
      
      % Test for exit.
      optimal = stats.pgNorms(end) < pgTol;
      smallf  = fy < fTol;
      if optimal || smallf
         break;
      end
      
   end
   x = stats.xOut;
   stats.iterations = k;
   stats.time = toc(stats.time);
end

function [y,fy,gy,ny,ax,opts] = linesearch(FG,P,x,fx,gx,nx,ax,opts)
   rdot = @(X,Y)real(dot(X(:),Y(:)));
   if ~isfield(opts,'armijo') || isempty(opts.armijo)
      opts.armijo = 1e-4;
   end
   if ~isfield(opts,'wolfe') || isempty(opts.wolfe)
      opts.wolfe = 0.9;
   end
   if ~isfield(opts,'rhodecrease') || isempty(opts.rhodecrease)
      opts.rhodecrease = (sqrt(5)+1)/2;
   end
   if ~isfield(opts,'rhoincrease') || isempty(opts.rhoincrease)
      opts.rhoincrease = (sqrt(5)+1)/2;
   end
   if ~isfield(opts,'minstep') || isempty(opts.minstep)
      opts.minstep = 2^-20;
   end
   if ~isfield(opts,'maxstep') || isempty(opts.maxstep)
      opts.maxstep = 2^20;
   end
   if ~isfield(opts,'maxit') || isempty(opts.maxit)
      opts.maxit = 2^4;
   end
   if ~isfield(opts,'F') || isempty(opts.F)
      opts.F = fx;
   end
   if ~isfield(opts,'kappamin') || isempty(opts.kappamin)
      opts.kappamin = 0;
   end
   if ~isfield(opts,'kappamax') || isempty(opts.kappamax)
      opts.kappamax = 1;
   end
   if ~isfield(opts,'kappa') || isempty(opts.kappa)
%       opts.kappa = 0.85;
      opts.kappa = 1;
   end
   if ~isfield(opts,'K') || isempty(opts.K)
      opts.K = 1;
   end
   y = P(x-ax*gx); [fy,gy] = FG(y); ny = (x-y)/ax-gx;
   exy   = fy-fx-rdot(gx,y-x);
   cxy   = rdot(gy-gx,y-x);
   gxymx = rdot(gx,y-x);
%    exy   = fy-fx-rdot(gx+nx,y-x);
%    cxy   = rdot(gy+ny-gx-nx,y-x);
%    gxymx = rdot(gx+nx,y-x);
   emax  = opts.F + (opts.armijo-1) * gxymx;
   cmin  = (opts.wolfe-1) * gxymx;
   if exy > emax
      decrease  = true;
      endsearch = (ax <= (opts.rhodecrease * opts.minstep));
   elseif cxy < cmin
      decrease  = false;
      endsearch = ((ax * opts.rhoincrease) >= opts.maxstep);
   else
      endsearch = true;
   end
   k = 1;
   opts.functions = fy;
   while ~endsearch && (k < opts.maxit)
      if decrease
         ax   = ax/opts.rhodecrease; y = P(x-ax*gx); [fy,gy] = FG(y); ny = (x-y)/ax-gx;
         exy  = fy-fx-rdot(gx,y-x);
         emax = opts.F + (opts.armijo-1)*rdot(gx,y-x);
%          exy  = fy-fx-rdot(gx+nx,y-x);
%          emax = opts.F + (opts.armijo-1)*rdot(gx+nx,y-x);
         endsearch = (ax < (opts.rhodecrease * opts.minstep)) || (exy <= emax);
      else
         axT   = ax*opts.rhoincrease; yT = P(x-axT*gx); [fyT,gyT] = FG(yT); nyT = (x-yT)/axT-gx;
         exy   = fyT-fx-rdot(gx,yT-x);
         cxy   = rdot(gyT-gx,yT-x);
         gxymx = rdot(gx,yT-x);
%          exy   = fyT-fx-rdot(gx+nx,yT-x);
%          cxy   = rdot(gyT-nyT-gx-nx,yT-x);
%          gxymx = rdot(gx+nx,yT-x);
         emax  = opts.F+(opts.armijo-1)*gxymx;
         cmin  = (opts.wolfe-1)*gxymx;
         endsearch = ((axT * opts.rhoincrease) > opts.maxstep) || (exy > emax) || (cxy >= cmin);
         if exy <= emax
            ax = axT; y = yT; fy = fyT; gy = gyT; ny = nyT;
         end
      end
      k = k + 1;
      opts.functions = [opts.functions;fy];
   end
   opts.K = 1+opts.kappa*opts.K;
   opts.F = opts.F+(fy-opts.F)/opts.K;
   opts.iterations = k;
end

% function [y,fy,gy,ny,ax,opts] = linesearch(FG,P,x,fx,gx,nx,ax,opts)
%    rdot = @(X,Y)real(dot(X(:),Y(:)));
%    if ~isfield(opts,'armijo') || isempty(opts.armijo)
%       opts.armijo = 1e-4;
%    end
%    if ~isfield(opts,'wolfe') || isempty(opts.wolfe)
%       opts.wolfe = 0.9;
%    end
%    if ~isfield(opts,'rhodecrease') || isempty(opts.rhodecrease)
%       opts.rhodecrease = (sqrt(5)+1)/2;
%    end
%    if ~isfield(opts,'rhoincrease') || isempty(opts.rhoincrease)
%       opts.rhoincrease = (sqrt(5)+1)/2;
%    end
%    if ~isfield(opts,'minstep') || isempty(opts.minstep)
%       opts.minstep = 2^-20;
%    end
%    if ~isfield(opts,'maxstep') || isempty(opts.maxstep)
%       opts.maxstep = 2^20;
%    end
%    if ~isfield(opts,'maxit') || isempty(opts.maxit)
%       opts.maxit = 2^4;
%    end
%    if ~isfield(opts,'F') || isempty(opts.F)
%       opts.F = fx;
%    end
%    if ~isfield(opts,'kappamin') || isempty(opts.kappamin)
%       opts.kappamin = 0;
%    end
%    if ~isfield(opts,'kappamax') || isempty(opts.kappamax)
%       opts.kappamax = 1;
%    end
%    if ~isfield(opts,'kappa') || isempty(opts.kappa)
%       opts.kappa = .85;
%    end
%    if ~isfield(opts,'K') || isempty(opts.K)
%       opts.K = 1;
%    end
%    y = P(x-ax*gx); [fy,gy] = FG(y); ny = (x-y)/ax-gx;
%    exy   = fy-fx-rdot(gx,y-x);
%    cxy   = rdot(gy-gx,y-x);
%    gxymx = rdot(gx,y-x);
%    emax  = opts.F + (opts.armijo-1) * gxymx;
%    cmin  = (opts.wolfe-1) * gxymx;
%    if exy > emax
%       decrease  = true;
%       endsearch = (ax <= (opts.rhodecrease * opts.minstep));
%    elseif cxy < cmin
%       decrease  = false;
%       endsearch = ((ax * opts.rhoincrease) >= opts.maxstep);
%    else
%       endsearch = true;
%    end
%    k = 1;
%    opts.functions = fy;
%    dx = y-x; ax = 1;
%    while ~endsearch && (k < opts.maxit)
%       if decrease
%          ax   = ax/opts.rhodecrease; y = x+ax*dx; [fy,gy] = FG(y); ny = (x-y)/ax-gx;
%          exy  = fy-fx-rdot(gx,y-x);
%          emax = opts.F + (opts.armijo-1)*rdot(gx,y-x);
%          endsearch = (ax < (opts.rhodecrease * opts.minstep)) || (exy <= emax);
%       else
%          axT   = ax*opts.rhoincrease; yT = x+axT*dx; [fyT,gyT] = FG(yT); nyT = (x-yT)/axT-gx;
%          exy   = fyT-fx-rdot(gx,yT-x);
%          cxy   = rdot(gyT-gx,yT-x);
%          gxymx = rdot(gx,yT-x);
%          emax  = opts.F+(opts.armijo-1)*gxymx;
%          cmin  = (opts.wolfe-1)*gxymx;
%          endsearch = ((axT * opts.rhoincrease) > opts.maxstep) || (exy > emax) || (cxy >= cmin);
%          if exy <= emax
%             ax = axT; y = yT; fy = fyT; gy = gyT; ny = nyT;
%          end
%       end
%       if isnan(fy) || ~isfinite(fy)
%          endsearch = true;
%       end
%       k = k + 1;
%       opts.functions = [opts.functions;fy];
%    end
%    opts.K = 1+opts.kappa*opts.K;
%    opts.F = opts.F+(fy-opts.F)/opts.K;
%    opts.iterations = k;
% end

function a = estimatestep(s,y,opts)
   if nargin == 2 || ~isstruct(opts)
      opts = [];
   end
   if ~isfield(opts,'estimator') || isempty(opts.estimator)
      opts.estimator = 'bb1';
   end
   switch lower(opts.estimator)
      case 'bb1'
         a = bb1step(s,y,opts);
      case 'bb2'
         a = bb2step(s,y,opts);
      otherwise
         error('[ERROR] Steplength estimator ''%s'' not supported!',opts.estimator);
   end
end

function a = bb1step(s,y,opts)
   ns = norm(s(:));
   ra = real(dot(y(:)/ns,s(:)/ns));
   if 1 < (ra * opts.minstep)
      a = opts.minstep;
   elseif 1 > (ra * opts.maxstep)
      a = opts.maxstep;
   else
      a = 1/ra;
   end
end

function a = bb2step(s,y,opts)
   ny = norm(y(:));
   a  = real(dot(y(:)/ny,s(:)/ny));
   if isnan(a) || a < 0 || a > opts.maxstep
      a = opts.maxstep;
   elseif a < opts.minstep
      a = opts.minstep;
   end
end
