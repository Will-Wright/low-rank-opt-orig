function addToNFFT2(n)
   global nfft2;
   if isempty(nfft2)
      nfft2 = n;
   else
      nfft2 = nfft2 + n;
   end
end
