function pattern = runiform(h,w,d)
   pattern = randn(h,w,d);
   pattern = bsxfun(@rdivide,pattern,realsqrt(sum(dot(pattern,pattern,1),2)/(h*w)));
end
