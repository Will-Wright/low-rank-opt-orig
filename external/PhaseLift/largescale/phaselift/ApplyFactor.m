function y = ApplyFactor(v,x,z,theta)
% - 
%
% Input:  
%

% Author:     Thomas Strohmer 
%


y = (1-theta)*(x*(x'*v)) + theta*(z*(z'*v));






% end of file
