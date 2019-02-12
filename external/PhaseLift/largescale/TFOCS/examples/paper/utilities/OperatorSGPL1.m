function x = OperatorSGPL1(z,md,A,At)
% wrapper for functions for SPGL
% Example usage: (from Jerome Bobin)
%   G = @(z,md) OperatorSGPL1(z,md,R,Rt);
% then...
%   [x_spgl,rspgl,gspgl,infosspl] = spgl1_test(x_nepact,G, b,[], delta, Rt(b));
%
% TFOCS version 1.0, code by Michael Grant (mcg@cvxr.com) and Stephen Becker (srbecker@caltech.edu)

if (md == 1) 
    x = A(z);
end
if (md == 2) 
    x = At(z);
end

