function [ op, inp_dims, otp_dims ] = linop_stack( linearF, inp_dims, otp_dims )

%LINOP_STACK    Stacked linear operators.
%    OP = LINOP_STACK( linearF ), where linearF is a cell vector or cell
%    matrix, returns a function handle for a linear operator that accepts
%    TFOCS_TUPLE objects as input or output, as appropriate, and applies
%    the various linear operators in block matrix fashion.
%
%    If linearF has more than one row, then the output in its forward mode
%    or its input in adjoint mode is a TFOCS_TUPLE object. If linearF has
%    more than one column, then the output in its adjoint mode or its input
%    in forward mode is a TFOCS_TUPLE object.

if ~isa( linearF, 'cell' ),
    error( 'First argument must be a cell array.' );
end
[ m, n ] = size( linearF );
if nargin < 2 || isempty( inp_dims ),
    inp_dims = cell( 1, n );
end
if nargin < 3 || isempty( otp_dims ),
    otp_dims = cell( 1, m );
end
rescan = zeros(2,0);
for i = 1 : m,
    otp_d = otp_dims{i};
    for j = 1 : n,
        inp_d = inp_dims{j};
        lF = linearF{i,j};
        sZ = [];
        if isempty(lF),
        elseif isa( lF, 'function_handle' ),
            sZ = lF([],0);
        elseif ~isnumeric( lF ),
            error( 'Entries should be real matrices or linear operators.'  );
        elseif ~isreal(lF),  % Why? we now handle A: C --> C
            error( 'Matrix entries must be real.' );
        elseif numel(lF) > 1,
            sZ = size(lF);
            linearF{i,j} = linop_matrix( lF );
        elseif lF == 0,
            linearF{i,j} = [];
        else
            if lF == 1,
                linearF{i,j} = @(x,mode)x;
            else
                linearF{i,j} = @(x,mode)lF*x;
            end
            if ~isempty(otp_d),
                sZ = { otp_d, otp_d };
            elseif ~isempty(inp_d),
                sZ = { inp_d, inp_d };
            else
                rescan(:,end+1) = [i;j];
            end
        end
        if isempty( sZ ),
            continue;
        elseif isnumeric( sZ ),
            sZ = { [sZ(2),1], [sZ(1),1] };
        end
        if isempty(inp_d),
            inp_d = sZ{1};
        elseif ~isequal(inp_d,sZ{1}) && ~isempty( sZ{1} ) % adding Oct 12
            error( 'Incompatible dimensions in element (%d,%d) of the linear operator matrix', i, j );
        end
        if isempty(otp_d),
            otp_d = sZ{2};
        elseif ~isequal(otp_d,sZ{2}),
            error( 'Incompatible dimensions in element (%d,%d) of the linear operator matrix', i, j );
        end
        inp_dims{j} = inp_d;
    end
    otp_dims{i} = otp_d;
end

%
% In some cases, we cannot resolve the dimensions on the first pass:
% specifically, those entries that represent scalar scaling operations.
% In those cases, we know that the input and output dimensions must be the
% same, but we may not have yet determined either in the first pass. So
% we rescan those entries until all ambiguities are resolved or until no
% further progress is made.
%

while ~isempty(rescan),
    rescan_o = rescan;
    rescan = zeros(2,0);
    for ij = rescan,
        i = ij(1); j = ij(2);
        lF = linearF{i,j};
        if isnumeric(lF) && numel(lF) == 1,
            if isempty(inp_dims{j}),
                if isempty(otp_dims{i}),
                    rescan(:,end+1) = [i;j];
                    continue;
                else
                    inp_dims{j} = otp_dims{i};
                end
            elseif isempty(otp_dims{i}),
                otp_dims{i} = inp_dims{j};
            elseif ~isequal( inp_dims{i}, otp_dims{j} ),
                error( 'Incompatible dimensions in element (%d,%d) of the linear operator matrix', i, j );
            end
        end
    end
    % Prevent infinite loops
    if numel(rescan) == numel(rescan_o),
        break;
    end
end

if m == 1 && n == 1,
    op = linearF{1,1};
    inp_dims = inp_dims{1};
    otp_dims = otp_dims{1};
    if isempty(op),
        op = @linop_identity;
    end
elseif m == 1,
    otp_dims = otp_dims{1};
    op = @(x,mode)linop_stack_row( linearF,  n,    { inp_dims, otp_dims }, x, mode );
elseif n == 1,
    inp_dims = inp_dims{1};
    op = @(x,mode)linop_stack_col( linearF,  m,    { inp_dims, otp_dims }, x, mode );
else
    op = @(x,mode)linop_stack_mat( linearF, [m,n], { inp_dims, otp_dims }, x, mode );
end

function y = linop_stack_row( linearF, N, dims, x, mode )
switch mode,
    case 0,
        y = dims;
    case 1,
        y = 0;
        x = cell( x );
        for j = 1 : N,
            lF = linearF{j};
            if ~isempty(lF), y = y + lF(x{j},1); end
        end
    case 2,
        y = cell(1,N);
        for j = 1 : N,
            lF = linearF{j};
            if ~isempty(lF), y{j} = lF(x,2); else y{j} = 0*x; end
        end
        y = tfocs_tuple( y );
end
        
function y = linop_stack_col( linearF, N, dims, x, mode )
switch mode,
    case 0,
        y = dims;
    case 1,
        y = cell(1,N);
        for j = 1 : N,
            lF = linearF{j};
            if ~isempty(lF), y{j} = lF(x,1); else y{j} = 0*x; end
        end
        y = tfocs_tuple( y );
    case 2,
        y = 0;
        x = cell( x );
        for j = 1 : N,
            lF = linearF{j};
            if ~isempty(lF), y = y + lF(x{j},2); end
        end
end

function y = linop_stack_mat( linearF, sZ, dims, x, mode )
switch mode,
    case 0,
        y = dims;
    case 1,
        x = cell( x );
        y = cell( 1, sZ(1) );
        for i = 1 : sZ(1),
            ans = 0;
            for j = 1 : sZ(2),
                lF = linearF{i,j};
                if ~isempty(lF), ans = ans + lF(x{j},1); end
            end
            y{i} = ans;
        end
        y = tfocs_tuple( y );
    case 2,
        x = cell( x );
        y = cell( 1, sZ(2) );
        for j = 1 : sZ(2),
            ans = 0;
            for i = 1 : sZ(1),
                lF = linearF{i,j};
                if ~isempty(lF), ans = ans + lF(x{i},2); end
            end
            y{j} = ans;
        end
        y = tfocs_tuple( y );
end

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
