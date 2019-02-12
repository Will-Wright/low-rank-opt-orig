function pattern = rademacher(h,w,d)
   pattern = randi(2,[h,w,d])*2-3;
end
