function Pattern = MakePattern1d(n,n_illumin,patterntype,opts)
% - 
%
% Usage:      
%
% Input:  
%             
% Output:    
%

switch lower(patterntype)
   case {'gaussian','gauss'}
      Pattern = randn(n,n_illumin) + i*randn(n,n_illumin);
      Pattern(:,1) = 1;
   case {'exponential','exp'}
      Pattern = exp(2*pi*i*rand(n,n_illumin));
      Pattern(:,1) = 1;
   case 'uniform'
      Pattern = rand(n,n_illumin);
      Pattern(:,1) = 1;
   case 'yonina'
      s = 1;
      t = [0:n-1]';
      h = ones(n,1); 
      Pattern = [h,1+exp(1i*2*pi*t*s/n),1-1i*exp(1i*2*pi*t*s/n)];
   case 'yoninagauss'
      s = 1;
      t = [0:n-1]';
      h = ones(n,1); 
      g = randn(n,1) + i*randn(n,1);
      Pattern = [g.*h,g.*(1+exp(1i*2*pi*t*s/n)),g.*(1-1i*exp(1i*2*pi*t*s/n))];
   case {'binary','bernoulli','0/1'}
      Pattern = round(rand(n,n_illumin));
      Pattern(:,1) = 1;
   case {'sine','grating'}
      Pattern = SinePattern(n,n_illumin);
   case 'ptycho'
      Pattern = PtychoPattern(n,n_illumin,8);
end

Pattern = normmat(Pattern)*sqrt(size(Pattern,1));






% end of file
