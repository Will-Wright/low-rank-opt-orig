%SOLVER_INITIALIZE	TFOCS helper script
%	Performs the initializations common to all of the first-order solvers

% Process the options structure and string/value pairs. First, we replace the
% default values with any values specified by the user in the 'opts' structure.
% We ignore any fields in the 'opts' structure that do not match ours, so that
% users can re-use the opts structure for other purposes.

odef = struct( ...
    'maxIts',     Inf   ,...
    'maxCounts',  Inf   ,...
    'countOps',   false ,... % adjusted if maxCounts < Inf
    'saveHist',   true  ,...
    'adjoint',    false ,...
    'saddle',     false ,...
    'tol',        1e-8  ,...
    'errFcn',     {{}}  ,...
    'stopFcn',    {{}}  ,...
    'printEvery', 100   ,...
    'maxmin',     1     ,...
    'beta',       0.5   ,...
    'alpha',      0.9   ,...
    'L0',         Inf   ,... % See below for the "true" default
    'Lexact',     Inf   ,...
    'mu',         0     ,...
    'fid',        1     ,...
    'stopCrit',   1     ,... % the TFOCS paper used stopCrit = 2
    'restart',    Inf   );

% Calling the solver with a no arguments returns the default options structure
if nargin < 1,
    opts = odef;
    out = [];
    x = opts;
    if nargout == 0,
        disp('====== TFOCS default options =======');
        disp( opts );
    end
    return % from this script only; the calling function does not return
end

% [smoothF, affineF, projectorF, x0, opts ] = deal(varargin{:});
% Check for incorrect types
F_types = {'function_handle','cell','double'};
assert( ismember(class(smoothF),F_types),'smoothF is of wrong type' );
    

% Process the options structure, merging it with the default options
def_fields = fieldnames(odef)';
if nargin > 4 && ~isempty(opts),
    use_fields = zeros(size(def_fields));
    opt_fields = fieldnames(opts)';
    for k = opt_fields,
        k = k{1};
        ndx = find(strcmpi(def_fields,k));
        if ~isempty(ndx)
            if ~isempty(opts.(k)) && ~use_fields(ndx),
                odef.(def_fields{ndx}) = opts.(k);
            end
            use_fields(ndx) = use_fields(ndx) + 1;
        else
            % Warn if the field is not found in our options structure
            warnState = warning('query','backtrace');
            warning off backtrace
            warning('TFOCS:OptionsStructure',' Found extraneous field "%s" in options structure', k );
            warning(warnState);
        end
    end
    % Warn if fields appear twice due to capitalization; e.g., maxIts/maxits
    if any(use_fields>1),
        warnState = warning('query','backtrace');
        warning off backtrace
        warning('TFOCS:OptionsStructure',' Some fieldnames of the options structure are identical up to capitalization: unpredictable behavior');
        warning(warnState);
    end
end

% The default value of L0 is set to Lexact, if it is supplied
if isinf(odef.L0),
    if isinf(odef.Lexact),
        if odef.beta >= 1,
            error( 'For a fixed step size, L0 or Lexact must be supplied.' );
        end
        odef.L0 = 1;
    else
        odef.L0 = odef.Lexact;
    end
end
% If maxCounts is set, set countOps to true
if any(odef.maxCounts<Inf),
    odef.countOps = true;
end

% Now move the options into the current workspace
for k = def_fields,
    assignin('caller',k{1},odef.(k{1}));
end
opts = odef;

%
% Smooth function, pass 1
%

if isempty( smoothF ),
	error( 'Must supply a smooth function specification.' );
elseif isa( smoothF, 'function_handle' ),
    smoothF = { smoothF };
elseif ~isa( smoothF, 'cell' ),
    error( 'smoothF must be a function handle, or a cell array of them.' );
end
n_smooth = numel(smoothF);
saddle_ndxs = 1 : n_smooth;

%
% Projector, pass 1
%

if isa( projectorF, 'function_handle' ),
    projectorF = { projectorF };
elseif ~isa( projectorF, 'cell' ) && ~isempty( projectorF ),
    error( 'projectorF must be a function handle, or a cell array of them.' );
end
n_proj = numel(projectorF);

%
% Linear functions and affine offsets
%

% If the affine operator is anything *but* a cell array, convert it to one.
maxmin = sign(maxmin);
if isempty( affineF )
    affineF = { @linop_identity };
elseif isnumeric( affineF ),
    if ndims(affineF) > 1,
        error( 'Multidimensional arrays are not permitted.' );
    end
    if numel(affineF) == 1,
        identity_linop = affineF == 1;
        affineF = { linop_scale( affineF ) };
    else
        identity_linop = false;
        affineF = { linop_matrix( affineF ) };
    end
elseif isa( affineF, 'function_handle' ),
    affineF = { affineF };
elseif ~isa( affineF, 'cell' ),
    error( 'Invalid affine operator specification.' );
end

% If adjoint mode is specified, temporarily transpose the affineF cell
% array so that the rows and columns match the number of smooth functions
% and projectors, respectively. Then verify size compatibility.
if adjoint, affineF = affineF'; end
[ m_aff, n_aff ] = size( affineF );
if all( m_aff ~= n_smooth + [0,1] ) || n_proj && all( n_aff ~= n_proj + [0,1] ),
    error( 'The affine operation matrix has incompatible dimensions.' );
elseif n_proj == 0,
    n_proj = max( 1, m_aff - 1 );
end
inp_dims = cell( 1, n_proj );
otp_dims = cell( 1, n_smooth );

% If an additional affine portion of the objective has been specified,
% create an additional smooth function to contain it.
if m_aff == n_smooth + 1,
    otp_dims{end+1} = [1,1];
    smoothF{end+1} = smooth_linear( maxmin ),pause
    n_smooth = n_smooth + 1;
    for k = 1 : n_proj,
        offX = affineF{end,k};
        if isempty(offX),
            affineF{end,k} = 0;
        elseif isa(offX,'function_handle'),
            if adjoint, pos = 'row'; else pos = 'column'; end
            error( 'The elements in the last %s must be constants.', pos );
        elseif isnumeric(offX) && numel(offX) == 1 && offX == 0,
            affineF{end,k} = 0;
        elseif nnz(offX),
            affineF{end,k} = linop_dot( affineF{end,k}, adjoint );
        end
    end
end

% If an affine offset has been specified, integrate those offsets into
% each smooth function, then remove that portion of the array.
if n_aff == n_proj + 1,
    for k = 1 : n_smooth,
        offX = affineF{k,end};
        if isempty(offX),
            continue;
        elseif isa(offX,'function_handle'),
            if adjoint, pos = 'column'; else pos = 'row'; end
            error( 'The elements in the last %s must be constant matrices.', pos );
        elseif isnumeric(offX) && numel(offX) == 1 && offX == 0,
            continue;
        else
            otp_dims{k} = size(offX);
            if nnz(offX),
                smoothF{k} = @(x)smoothF{k}( x + offX );
            end
        end
    end
    n_aff = n_aff - 1;
    affineF(:,end) = [];
end

% Transpose back, if necessary; then check dimensions
if adjoint,
    affineF = affineF';
    [ linearF, otp_dims, inp_dims ] = linop_stack( affineF, otp_dims, inp_dims );
    linearF = linop_adjoint( linearF );
else
    [ linearF, inp_dims, otp_dims ] = linop_stack( affineF );
end
identity_linop = isequal( linearF, @linop_identity );
square_linop = identity_linop || isequal( inp_dims, otp_dims );
adj_arr = [0,2,1];
if countOps,
    apply_linear = @(x,mode)solver_apply( 3, linearF, x, mode );
else
    apply_linear = linearF;
end

%
% Smooth function, pass 2: integrate scale, counting
%

smoothF = smooth_stack( smoothF );
if maxmin < 0,
    smoothF = tfunc_scale( smoothF, -1 );
end
if countOps,
    apply_smooth = @(x)solver_apply( 1:(1+(nargout>1)), smoothF, x );
else
    apply_smooth = smoothF;
end

%
% Projector, pass 2: supply default, stack it up, etc.
%

if isempty( projectorF ),
    n_proj = 0;
    projectorF = proj_Rn;
else
    projectorF = prox_stack( projectorF );
end
if countOps,
    apply_projector = @(varargin)solver_apply( 4:(4+(nargout>1)), projectorF, varargin{:} );
else
    apply_projector = projectorF;
end

%
% Initialize the op counts
%

if countOps,
    global tfocs_count___
    tfocs_count___ = [0,0,0,0,0];
    maxCounts = maxCounts(:)';
end

%
% Construct the common initial values
%

L = L0;
theta = Inf;
f_v_old = Inf;
x = []; A_x = []; f_x = Inf; C_x = Inf; g_x = []; g_Ax = [];
restart_iter = 0;
warning_lipschitz = 0;
backtrack_simple = true;
backtrack_tol = 1e-10;
backtrack_steps = 0;

%
% Try to determine the size of the input, and construct the initial point.
%

% Attempt 1: From x0 itself
zero_x0 = true;
if ~isempty( x0 ),
    if isa( x0, 'tfocs_tuple' )
        x0  = cell(x0);
    end
    if isa( x0, 'cell' ),
        n_x0 = numel( x0 );
        if n_x0 == 1,
            x0 = x0{1};
        else
            x0 = tfocs_tuple( x0 );
        end
    else
        n_x0 = 1;
    end
    if n_proj && n_proj ~= n_x0,
        error( 'Size conflict detected between the projector and x0.' );
    end
    zero_x0 = ~nnz(x0);
% Attempt 2: From the linear operator dimensions
elseif ~size_ambig( inp_dims ),
    x0 = tfocs_zeros( inp_dims );
elseif ~size_ambig( otp_dims ),
    A_x = tfocs_zeros( otp_dims );
    x0 = apply_linear( A_x, 2 );
end
if isempty( x0 ),
    error( 'Could not determine the dimensions of the problem. Please supply an explicit value for x0.' );
end
x = x0;
C_x = 0;
%if isinf( C_x ),
%    C_x = proj_tsvd( ApplyZhandle, dimx );
%    if isinf( C_x ),
%        zero_x0 = false;
%        [ C_x, x ] = proj_tsvd(ApplyZhandle, dimx );
%    end
%end
if isempty( A_x ),
    if identity_linop || zero_x0 && square_linop,
        A_x = x;
    elseif ~zero_x0 || size_ambig( otp_dims ),
        A_x = apply_linear( x, 1 );
    else
        A_x = tfocs_zeros( otp_dims );
    end
end
% Final size check
[ isOK1, inp_dims ] = size_compat( size(x0), inp_dims );
[ isOK2, otp_dims ] = size_compat( size(A_x), otp_dims );
if ~isOK1 || ~isOK2,
    error( 'Could not determine the dimensions of the problem. Please supply an explicit value for x0.' );
end
[ f_x, g_Ax ] = smoothF( A_x );
if isinf( f_x ),
    error( 'The initial point lies outside of the domain of the smooth function.' );
end

% Theta advancement function
if mu > 0 && Lexact > mu,
    theta_scale = sqrt(mu / Lexact);
    theta_scale = ( 1 - theta_scale ) / ( 1 + theta_scale );
    advance_theta = @(theta_old,L,L_old) min(1,theta_old*theta_scale);
else
    advance_theta = @(theta_old,L,L_old) 2/(1+sqrt(1+4*(L/L_old)/theta_old.^2));
end

% Preallocate the arrays for the output structure
out.alg = alg; 
out.algorithm = algorithm;
if ~isempty(errFcn) && ~iscell(errFcn),
   errFcn = { errFcn };
end
errs = zeros(1,length(errFcn));
if nargout == 1,
    saveHist = false;
end
if saveHist,
    [ out.f, out.normGrad, out.stepsize, out.theta ] = deal( zeros(0,1) );
    if countOps,
        out.counts = zeros(0,length(tfocs_count___));
    end
    if ~isempty(errFcn),
        out.err = zeros(0,length(errs));
    end
end
if saddle,
    out.dual = [];
end
n_iter = 0;
status = '';

% Initialize the iterate values
y    = x;    z    = x;
A_y  = A_x;  A_z  = A_x;
C_y  = C_x;  C_z  = C_x;
f_y  = f_x;  f_z  = f_x;
g_y  = g_x;  g_z  = g_x;
g_Ay = g_Ax; g_Az = g_Ax;
norm_x = sqrt( tfocs_normsq( x ) );

% Special setup for constant step sizes
if beta >= 1,
    beta = 1;
	alpha = 1;
end

% Print the opening text
if fid && printEvery,
	fprintf(fid,'%s\n',algorithm);
	fprintf(fid,'Iter    Objective   |dx|/|x|    step');
    if countOps, fprintf( fid, '       F     G     A     N     P' ); end
    if ~isempty(errFcn), fprintf( fid, '      errors' ); end
    fprintf( fid, '\n' );
	fprintf(fid,'----+----------------------------------' );
    if countOps, fprintf( fid, '+-------------------------------' ); end
    if ~isempty(errFcn), fprintf( fid, '+-------------------' ); end
    fprintf(fid,'\n');
end

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.

