function addToNIFFT2(n)
   global nifft2;
   if isempty(nifft2)
      nifft2 = n;
   else
      nifft2 = nifft2 + n;
   end
end
