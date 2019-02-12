function y = FactoredIllumination(X,Illumin)
%  FactoredIllumination computes data as acquired by diffracted structured
%  illumination patterns, but X is already assumed to be in factorized
%  form so that X*X' = Z, where Z is an hermitian pos.semidef. low-rank matrix
%  with [V,D] = eig(Z) and X = V*D^(1/2).
%
% Usage:      y = FactoredIllumination(X,Illumin)
%

% Author:     Thomas Strohmer 
%

% Track number of function calls.
persistent nCalls
if ~exist('nCalls','var')
   nCalls = 0;
end
if nargin == 0
   y = nCalls;
   nCalls = 0;
   return
end
nCalls = nCalls + 1;


if length(size(Illumin)) == 2;
  [n,L] = size(Illumin);

  [~,RNK] = size(X);

  y = zeros(n,L);
  for k = 1:L;
      D = Illumin(:,k);
      for j = 1:RNK
          y(:,k) = y(:,k) + abs(fft(D(:).*X(:,j))).^2;
      end
  end
  y = y/n;

elseif length(size(Illumin)) == 3;

  [m1,m2,m3] = size(Illumin);

  [n,RNK] = size(X);
  n1 = round(sqrt(n));
  n2 = n1;

  for k = 1:m3;
      M = Illumin(:,:,k);
      ypart = zeros(m1,m2);
      for j = 1:RNK
         ypart = ypart + abs(fft2(reshape(M(:).*X(:,j),n1,n2))).^2;
      end
      y(:,k) = ypart(:)/(m1*m2);
  end
else 
  error('Wrong array dimensions');
end




% end of file
