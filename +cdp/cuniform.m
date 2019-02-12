function pattern = cuniform(h,w,d)
   pattern = complex(randn(h,w,d),randn(h,w,d));
   pattern = bsxfun(@rdivide,pattern,realsqrt(sum(dot(pattern,pattern,1),2)/(h*w)));
end
