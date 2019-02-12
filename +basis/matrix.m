classdef matrix
   properties
      support
      m1
      m2
      m
      n
   end
   methods
      function self = matrix(support)
         self.support = support;
         self.m1      = size(self.support,1);
         self.m2      = size(self.support,2);
         self.m       = self.m1*self.m2;
         self.n       = nnz(self.support);
      end
      function out = forward(self,x)
         out = zeros(self.m1,self.m2,numel(x)/self.n,'like',x);
         out(self.support(:,:,ones(size(out,3),1))) = x(:);
      end
      function out = adjoint(self,x)
         k   = numel(x)/self.m;
         x   = reshape(x,[self.m1,self.m2,k]);
         out = reshape(x(self.support(:,:,ones(k,1))),self.n,k);
      end
      function sz = size(self,dim)
         sz = [self.m,self.n];
         if nargin == 2
            if dim > 2
               sz = 1;
            elseif dim < 1
               error('matrix:size',['Error using <strong>size</strong>\n'  ...
                     'Dimension argument must be a positive integer\n' ...
                     'scalar within indexing range.']);
            else
               sz = sz(dim);
            end
         end
      end
   end
end