function pattern = cgaussian(h,w,d)
   pattern = complex(randn(h,w,d),randn(h,w,d))/sqrt(2);
end
