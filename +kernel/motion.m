function kernel = motion(sz,len,theta)
   K      = fspecial('motion',len,theta);
   sza    = 1+ceil ((sz-size(K)-1)/2);
   szb    =   floor((sz-size(K)-1)/2);
   kernel = full(fftshift(blkdiag(sparse(sza(1),sza(2)),K,sparse(szb(1),szb(2)))));
end
