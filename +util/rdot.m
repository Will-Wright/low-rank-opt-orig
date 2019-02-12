function rdotxy = rdot(x,y)
   % RDOT computes the real-part of a vectorized complex inner-product
   rdotxy = real(dot(x(:),y(:)));
end
