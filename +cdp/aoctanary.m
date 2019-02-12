function pattern = aoctanary(h,w,d)
   sz      = [h,w,d];
   b1      = [-1;1;-1i;1i];
   b2      = [1;1;1;1;realsqrt(6)];
   pattern = b1(randi(4,sz,'uint8')).*b2(randi(5,sz,'uint8'));
end
