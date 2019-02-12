function kernel = cmotion(sz,r)
   K      = zeros(2*ceil(r)+1);
   K      = double(bsxfun(@(kx,ky)double(abs(realsqrt(kx.^2+ky.^2)-norm(size(K))/2)<=1).*exp(-(realsqrt(kx.^2+ky.^2)-norm(size(K))/2).^2),(1:size(K,1))',1:size(K,2))); K = K/sum(K(:));
   sza    = 1+ceil ((sz-size(K)-1)/2);
   szb    =   floor((sz-size(K)-1)/2);
   kernel = full(fftshift(blkdiag(sparse(sza(1),sza(2)),K,sparse(szb(1),szb(2)))));
end
