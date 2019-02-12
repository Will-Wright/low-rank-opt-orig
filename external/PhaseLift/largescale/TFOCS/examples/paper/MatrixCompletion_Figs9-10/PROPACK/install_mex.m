function install_mex(VERB)
% install_mex
%   installs mex files for SVT and the PROPACK subpackage
% install_mexSimple( verbose )
%   will display more verbose diagnoistic messages
%   if the variable "verbose" is true.
%
% Stephen Becker, Jan 2010, SVT package. srbecker@caltech.edu

SVT = false; % do not install special SVT files

if nargin < 1 || ~VERB
    VERBOSE = '';
else
    VERBOSE = '-v';
end

X='';
if ispc
    cc = get_compiler_config();
    if strcmpi(cc,'microsoft'),         cc = 'microsoft';
    else,         cc = 'lcc';     % Matlab's bundled compiler
    end
    c = computer;
    if strfind(c,'64')
        libpath = fullfile(matlabroot,'extern','lib','win64',cc);
    else
        libpath = fullfile(matlabroot,'extern','lib','win32',cc);
    end
    LAPACK = fullfile(libpath,'libmwlapack.lib');
    BLAS = fullfile(libpath,'libmwblas.lib');
    WIN = '-DWINDOWS';
else
    WIN = '-UWINDOWS';
    LAPACK = '-lmwlapack';
    % on linux/unix, sometimes MATLAB won't install its mwblas library
    % if the system already has a blas library.  So check for this:
    % (note, not in the same location as on Windows)
    c = lower(computer);
    if ismac, suffix = '.dylib'; else suffix = '.so'; end
    blasFile = fullfile(matlabroot,'bin',c, ['libmwblas',suffix] );
    if exist(blasFile,'file')
        BLAS = '-lmwblas';
    else
        BLAS = '-lblas';
        X='-DNO_BLAS';  % tell it not to include blas.h
    end
end
    

EXT = [];  % for native.  Use this when compiling for your own computer
% 2006a is v 7.2
if verLessThan('matlab', '7.3')
    LARGEARRAYDIMS = [];
    Y='-DNO_MATRIX_H'; % don't include matrix.h 'cause it doesn't exist!
else
    LARGEARRAYDIMS = '-largeArrayDims';
    Y=[];
end

OPT = '-O';

% compile PROPACK
mexHelper(VERBOSE,WIN,OPT,LARGEARRAYDIMS,EXT,'bdsqr_mex.c','dbdqr.c','-output','bdsqr',LAPACK,BLAS,X,Y);
mexHelper(VERBOSE,WIN,OPT,LARGEARRAYDIMS,EXT,'reorth_mex.c','reorth.c','-output','reorth',LAPACK,BLAS,X,Y);
mexHelper(VERBOSE,    OPT,LARGEARRAYDIMS,EXT,'computeInt.c');
% compile SVT code
if SVT
 mexHelper(VERBOSE,WIN,OPT,LARGEARRAYDIMS,EXT,'XonOmega.c',BLAS,X,Y);
 mexHelper(VERBOSE,WIN,OPT,LARGEARRAYDIMS,EXT,'XonOmegaTranspose.c',BLAS,X,Y);
 mexHelper(VERBOSE,WIN,OPT,LARGEARRAYDIMS,EXT,'updateSparse.c',BLAS,X,Y);
 mexHelper(VERBOSE,WIN,OPT,LARGEARRAYDIMS,EXT,'smvp.c',Y);
end

fprintf('If you see this, the installation was probably successful\n');
if SVT
 fprintf('Please run "test_MEX.m" and "test_PROPACK.m" to verify your installation\n\n');
else
 fprintf('Please run "test_PROPACK.m" to verify your installation\n\n');
end


function cc = get_compiler_config()
    % tested on Windows w/ R2008 only
    % This has to be in a function, otherwise old versions of matlab
    % get confused because "mex" is used as a structure (well, a class)
    % AND as a function.
    try 
        % this requires both a new verson of matlab and
        % that a compiler has been selected
        cc = mex.getCompilerConfigurations('C');
        cc = cc.Manufacturer;
        % Watch out for this error later (in old versions of matlab)
%         MATLAB:mir_error_function_previously_indexed_by_dot
    catch
        cc = [];
        fprintf('You may want to run ''mex -setup'' to setup the mex compiler,\n if you''ve never used the mex compiler before\n');
        wbsite='http://www.mathworks.com/support/solutions/en/data/1-6IJJ3L/index.html?solution=1-6IJJ3L';
        fprintf('If you have version 2008b or newer, and 64-bit Windows,\n');
        fprintf('then MATLAB does not come with a builtin compiler\n');
        fprintf('If you need a free C/C++ compiler, please see this mathworks website:\n%s\n',wbsite);
    end


function mexHelper(varargin)
  n = length(varargin);
  indx = [];
  for i = 1:n
      if ~isempty( varargin{i} )
          indx = [indx,i];
      end
  end
  mex( varargin{indx} )
