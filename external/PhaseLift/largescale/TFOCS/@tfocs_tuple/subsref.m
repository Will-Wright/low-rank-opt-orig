function varargout = subsref( x, varargin )
[ varargout{1:nargout} ] = subsref( x.value_, varargin{:} );

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
