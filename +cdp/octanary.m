function pattern = octanary(h,w,d)
   sz      = [h,w,d];
   b1      = [-1;1;-1i;1i];
   b2      = [repmat(sqrt(0.5),4,1); sqrt(3)];
   pattern = b1(randi(4,sz)) .* b2(randi(5,sz));
end
