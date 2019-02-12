function pattern = normalize(pattern)
   norms   = realsqrt(dot(pattern,pattern,3));
   zeros   = norms < eps;
   pattern = bsxfun(@rdivide,pattern,norms);
   if any(zeros(:))
      warning('[cdp.normalize] Pattern contains %d zero entries!',nnz(zeros));
      pattern(zeros) = 0;
   end
end
