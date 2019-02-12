% Implementation of the Wirtinger Flow (WF) algorithm presented in the paper 
% "Phase Retrieval via Wirtinger Flow: Theory and Algorithms" 
% by E. J. Candes, X. Li, and M. Soltanolkotabi

% The input data are coded diffraction patterns about an RGB image. Each
% color band is acquired separately (the same codes are used for all 3 bands)

% Code written by M. Soltanolkotabi and E. J. Candes

%%  Read Image 

% Below X is n1 x n2 x 3; i.e. we have three n1 x n2 images, one for each of the 3 color channels  
namestr = 'stanford' ;
stanstr = 'jpg'      ;
X       = mat2gray(imread([namestr,'.',stanstr])) ;
n1      = size(X,1)                               ;
n2      = size(X,2)                               ;

%% Make masks and linear sampling operators  

% Each mask has iid entries following the octanary pattern in which the entries are 
% distributed as b1 x b2 where 
% b1 is uniform over {1, -1, i, -i} (phase) 
% b2 is equal to 1/sqrt(2) with prob. 4/5 and sqrt(3) with prob. 1/5 (magnitude)

rng('default');

L = 21;                  % Number of masks  
Masks = zeros(n1,n2,L);  % Storage for L masks, each of dim n1 x n2

% Sample phases: each symbol in alphabet {1, -1, i , -i} has equal prob. 
for ll = 1:L, Masks(:,:,ll) = util.randsrc(n1,n2,[1i -1i 1 -1]); end

% Sample magnitudes and make masks 
temp = rand(size(Masks));
Masks = Masks .* ( (temp <= 0.2)*sqrt(3) + (temp > 0.2)/sqrt(2) );

% Make linear operators; 
A = @(I)  fft2(conj(Masks) .* reshape(repmat(I,[1 L]), size(I,1), size(I,2), L));  % Input is n1 x n2 image, output is n1 x n2 x L array
At = @(Y) sum(Masks .* ifft2(Y), 3) * size(Y,1) * size(Y,2);                       % Input is n1 x n2 x L array, output is n1 x n2 image

%% Prepare structure to save intermediate results 

ttimes   = [150,300];        % Iterations at which we will save info 
ntimes = length(ttimes)+1;   % +1 because we will save info after the initialization 
Xhats    = cell(1,ntimes);
for mm = 1:ntimes, Xhats{mm} = zeros(size(X));  end
Times    = zeros(3,ntimes);

%% Wirtinger flow  

npower_iter = 50;                   % Number of power iterations 
T = max(ttimes);                    % Max number of iterations
tau0 = 330;                         % Time constant for step size 
mu = @(t) min(1-exp(-t/tau0), 0.4); % Schedule for step size 

for rgb = 1:3, 
    fprintf('Color band %d\n', rgb)
    x = squeeze(X(:,:,rgb)); % Image x is n1 x n2 
    Y = abs(A(x)).^2;        % Measured data 
    
    % Initialization
    z0 = randn(n1,n2); z0 = z0/norm(z0,'fro'); % Initial guess 
    tic
    for tt = 1:npower_iter, 
        z0 = At(A(z0)); z0 = z0/norm(z0,'fro');
    end
    Times(rgb,1) = toc;
    
    normest = sqrt(sum(Y(:))/numel(Y)); % Estimate norm to scale eigenvector  
    z = normest * z0;                   % Apply scaling 
    Xhats{1}(:,:,rgb) = exp(-1i*angle(trace(x'*z))) * z; % Initial guess after global phase adjustment 
 
    % Loop    
    fprintf('Done with initialization, starting loop\n')
    tic
    for t = 1:T,
        Bz = A(z);
        C  = (abs(Bz).^2-Y) .* Bz;
        grad = At(C)/numel(C);                   % Wirtinger gradient            
        z   = z - mu(t)/normest^2 * grad;        % Gradient update 
        
        ind =  find(t == ttimes);                % Store results 
        if ~isempty(ind), 
             Xhats{ind+1}(:,:,rgb) = exp(-1i*angle(trace(x'*z))) * z; 
             Times(rgb,ind+1) = toc;
        end
        fprintf('%4i  %10.2e\n',t, norm(grad));
    end
    
end
fprintf('All done!\n')

%% Show some results 

iter = [0 ttimes];
Relerrs = zeros(1,ntimes);
for mm = 1:ntimes; 
    fprintf('Mean running times after %d iterations: %.1f\n', iter(mm), mean(Times(:,mm)))
    Relerrs(mm) = norm(Xhats{mm}(:)-X(:))/norm(X(:)); 
    fprintf('Relative error after %d iterations: %f\n', iter(mm), Relerrs(mm))  
    fprintf('\n')
end

for tt = 1:ntimes, 
    figure; imshow(mat2gray(abs(Xhats{tt})),[]);
end