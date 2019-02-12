classdef haar
   properties
      support
      m1
      m2
      m
      n
      h
   end
   methods
      function self = haar(support)
         self.support = support;
         self.m1      = size(self.support,1);
         self.m2      = size(self.support,2);
         self.m       = self.m1*self.m2;
         self.n       = nnz(self.support);
         self.h       = repmat(realsqrt(.5),1,2);
      end
      function out = forward(self,x)
         k   = numel(x)/self.n;
         x   = reshape(x,self.n,k);
         out = zeros(self.m1,self.m2,k,'like',x);
         for i = 1:k
            out(:,:,i) = self.atomicForward(x(:,i));
         end
      end
      function out = adjoint(self,x)
         k   = numel(x)/self.m;
         x   = reshape(x,[self.m1,self.m2,k]);
         out = zeros(self.n,k,'like',x);
         for i = 1:k
            out(:,i) = self.atomicAdjoint(x(:,:,i));
         end
      end
      function sz = size(self,dim)
         sz = [self.m,self.n];
         if nargin == 2
            if dim > 2
               sz = 1;
            elseif dim < 1
               error('haar:size',['Error using <strong>size</strong>\n'  ...
                     'Dimension argument must be a positive integer\n' ...
                     'scalar within indexing range.']);
            else
               sz = sz(dim);
            end
         end
      end
   end
   methods ( Access = private )
      function out = atomicForward(self,x)
         out = zeros(self.m1,self.m2,'like',x);
         out(self.support) = x(:);
         if isreal(out)
            out = rwt.midwt(out,self.h);
            util.addToNIDWT(1);
         else
            out = complex(rwt.midwt(real(out),self.h),rwt.midwt(imag(out),self.h));
            util.addToNIDWT(2);
         end
      end
      function out = atomicAdjoint(self,x)
         x = reshape(x,self.m1,self.m2);
         if isreal(x)
            out = rwt.mdwt(x,self.h);
            util.addToNDWT(1);
         else
            out = complex(rwt.mdwt(real(x),self.h),rwt.mdwt(imag(x),self.h));
            util.addToNDWT(2);
         end
         out = out(self.support);
      end
   end
end
