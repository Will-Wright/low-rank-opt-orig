function y = AdjGradient(g_Ay, v, Illumin)
%
% y = AdjGradient(g_Ay, v, Illumin) computes (A^* g_Ay)v,
% which is needed for computing Augmented Neg. Loglikelihood.
%
% y = AdjGradient() returns the number of calls to this routine, and resets
% the counter.

% Track number of adjoint calls.
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
   
   %   y = zeros(n,1);
   %
   %   for k = 1:L;
   %       D = Illumin(:,k);
   %       vpart = fft(D.*v);
   %       vpart = vpart.*g_Ay(:,k);
   %       vpart = conj(D).*ifft(vpart);
   %       y = y + vpart;
   %   end
   
   % MPF, Jan 13: This is a vectorized version of the above.
   M = Illumin;
   tmp = bsxfun(@times, M, v);
   tmp = fft(tmp);
   tmp = bsxfun(@times, reshape(g_Ay, n, L), tmp);
   tmp = ifft(tmp);
   tmp = bsxfun(@times, conj(M), tmp);
   y = sum(tmp, 2);

  
elseif length(size(Illumin)) == 3;

  [~,~,L] = size(Illumin);

  nn = length(v);
  n = round(sqrt(nn));

  y = zeros(nn,1);

  for k = 1:L;
    D = Illumin(:,:,k);
    vpart = fft2(reshape(D(:).*v,n,n));
    vpart = vpart(:).*g_Ay(:,k);
    vpart = conj(D).*ifft2(reshape(vpart,n,n));
    y = y + vpart(:);
  end

else
  error('Wrong array dimensions');
end
