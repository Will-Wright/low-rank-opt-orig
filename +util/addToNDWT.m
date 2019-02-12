function addToNDWT(n)
   global ndwt;
   if isempty(ndwt)
      ndwt = n;
   else
      ndwt = ndwt + n;
   end
end
