function pattern = quaternary(h,w,d)
   alphabet = [-1;1;-1i;1i];
   pattern  = alphabet(randi(4,[h,w,d]));
end
