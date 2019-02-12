function  Y=rotrc(X,row,col)
% ROTRC    - rotation of rows and columns of a matrix
%
% USAGE:  Y=rotrc(X,row,col);
% 
% Input : X    ...   matrix to be rotated
%         row  ...   rotate rows by row (integer)
%         col  ...   rotate columns by col (integer), 
%                    default: col=0
% Output: Y    ...   rotated version of X
%         
% see also FFTSHIFT

% Author : H.Steier and H.G.Feichtinger, 21.Nov.1989 
% Mods.  : H.G.Feichtinger, 18.Feb.1991  
%          M.Rauth,         08.Jan.1995  
%
% Copyright (c) NUHAG, Dept.Math., University of Vienna, AUSTRIA
%    Permission is granted to modify and re-distribute this 
%    code in any manner as long as this notice is preserved. 
%    All standard disclaimers apply. 
% 
if nargin<3
   col=0; 
end; 
[hig,wid]=size(X);
% 
% rotate rows
% 
r1=hig-rem(row,hig);
if r1>=hig;  
   r1=r1-hig; 
end;
Y=[X((r1+1):hig,:);X(1:r1,:)];
% 
% rotate columns
% 
r1=wid-rem(col,wid);
if r1>=wid;  
   r1=r1-wid; 
end;
Y=[Y(:,(r1+1):wid),Y(:,1:r1)];
