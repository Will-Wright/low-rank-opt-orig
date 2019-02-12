function Pattern = PtychoPattern(pattern,shiftpar);

%%% returns n x n_illumin matrix, each column being an illumination pattern 
%%% which modulates the signal

[n,m] = size(pattern);

if m > 1
  num_n = round(n/shiftpar(1));
  num_m = round(m/shiftpar(2));
  for k = 1:num_n 
    for j = 1:num_m
      ndx = (k-1)*num_m + j;
      Pattern(:,:,ndx) = rotrc(pattern,(k-1)*shiftpar(1),(j-1)*shiftpar(2));
    end
  end
else
  num_n = round(n/shiftpar);
  for k = 1:num_n 
    Pattern(:,k) = rot(pattern,(k-1)*shiftpar);;
  end
end

