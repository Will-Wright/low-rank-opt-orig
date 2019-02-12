function [z, info] = wflow(M, X0, varargin)
ip = inputParser;
ip.KeepUnmatched = true;
ip.addParameter('verbosity', true);
ip.addParameter('iterations', 300);
ip.addParameter('gTol', 1e-6);
ip.addParameter('Y', []);
ip.parse(varargin{:});
verbosity = ip.Results.verbosity;
iterations = ip.Results.iterations;
gTol = ip.Results.gTol;

% Implementation of the Wirtinger Flow (WF) algorithm presented in the paper 
% "Phase Retrieval via Wirtinger Flow: Theory and Algorithms" 
% by E. J. Candes, X. Li, and M. Soltanolkotabi

% The input data are coded diffraction patterns about an RGB image. Each
% color band is acquired separately (the same codes are used for all 3 bands)

% Code written by M. Soltanolkotabi and E. J. Candes
% 29 Jan 2015: encapsulate in function; MPF.


% Extract sizes from mask.
[n1, n2, L] = size(M);


% Make linear operators
A = @(I)  fft2(conj(M) .* reshape(repmat(I,[1 L]), size(I,1), size(I,2), L));  % Input is n1 x n2 image, output is n1 x n2 x L array
At = @(Y) mean(M .* ifft2(Y), 3) ;                       % Input is n1 x n2 x L array, output is n1 x n2 image

% Measure the data.
if isempty(ip.Results.Y)
   Y = abs(A(X0)).^2;
else
   Y = ip.Results.Y;
end

% Wirtinger flow  
npower_iter = 50;                   % Number of power iterations 
tau0 = 330;                         % Time constant for step size 
mu = @(t) min(1-exp(-t/tau0), 0.4); % Schedule for step size 
    
% Initialization
z0 = randn(n1,n2); z0 = z0/norm(z0,'fro'); % Initial guess
tStart = tic;
for tt = 1:npower_iter,
   z0 = At(Y.*A(z0)); z0 = z0/norm(z0,'fro');
end

normest = sqrt(sum(Y(:))/numel(Y)); % Estimate norm to scale eigenvector
z = normest * z0;                   % Apply scaling
 
% Loop
printf('Done with initialization, starting loop\n')
printf('%4s  %10s  %10s\n','iter','rNorm','gNorm');
t = 0;
while true
   t = t + 1;
   Bz = A(z);
   R = abs(Bz).^2-Y;
   C = R .* Bz;
   grad = At(C);            % Wirtinger gradient
   z = z - mu(t)/normest^2 * grad;   % Gradient update   
   gNorm = norm(grad(:));            % MPF edit: vec grad (not op norm)
   rNorm = norm(R(:));

   printf('%4i  %10.2e  %10.2e\n', t, rNorm, gNorm);
   if t == 1
      gTol0 = gTol*(1+gNorm);
   end
   

   if gNorm <= gTol0
      status = 'feasible';
      break
   elseif t >= iterations
      status = 'iterations';
      break
   end
   
end
    
printf('All done!\n');

info.time = toc(tStart);
info.status = status;
info.nAdjoint =  npower_iter + t;
info.nMeasure =  npower_iter + t;
info.nfft = (npower_iter + t)*L*2;
info.iterations = t;


% Local function (shared workspace)

   function printf(varargin)
      if verbosity
         fprintf(varargin{:});
      end
   end

end % function wflow
