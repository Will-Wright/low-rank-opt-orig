function addToNIDWT(n)
   global nidwt;
   if isempty(nidwt)
      nidwt = n;
   else
      nidwt = nidwt + n;
   end
end
