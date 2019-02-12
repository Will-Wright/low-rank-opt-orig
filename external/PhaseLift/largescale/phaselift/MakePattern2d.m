function Pattern = MakePattern2d(n,m,n_illumin,patterntype,opts)
% MakePattern2d - makes 2D masks for structured Illuminations
%
% Usage:   Pattern = MakePattern2d(n,m,n_illumin,patterntype,opts)
%
% Input:  
%             
% Output:    
%

switch lower(patterntype)
   case {'gaussian','gauss'}
      Pattern = randn(n,m,n_illumin) + i*randn(n,m,n_illumin);
   case {'exponential','exp'}
      Pattern = exp(2*pi*i*rand(n,m,n_illumin));
   case 'uniform'
      Pattern = rand(n,m,n_illumin);
   case {'binary','bernoulli','0/1'}
      Pattern = round(rand(n,m,n_illumin));
   case 'ptycho'
     Pattern = PtychoPattern(n,m,n_illumin,8);
end






% end of file
