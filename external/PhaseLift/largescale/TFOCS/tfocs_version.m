function tfocs_ver2 = tfocs_version

%TFOCS_VERSION   Version information.
%    Prints information about the version of TFOCS being used, and the
%    system on which it is running. When submitting a bug report, please
%    run this function and include its output in your report.

tfocs_ver = '1.0a';
if nargout == 0,
    fprintf( 'TFOCS version %s\n', tfocs_ver );
    verd = ver('MATLAB');
    fprintf( 'MATLAB version %s %s on %s\n', verd.Version, verd.Release, computer );
else
    tfocs_ver2 = tfocs_ver;
end

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.