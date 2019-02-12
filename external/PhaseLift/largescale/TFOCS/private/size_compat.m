function [ a, sX ] = size_compat( sX, sY )
a = true;
switch class( sX ),
    case 'double',
        if isempty( sX ) || all( sX == 0 ),
            sX = sY;
        elseif isempty( sY ) || all( sY == 0 ),
        elseif ~isequal( sX, sY ),
            a = false;
        end
    case 'cell',
        if ~isa( sY, 'cell' ) || numel( sX ) ~= numel( sY ) || isa( sX{1}, 'function_handle' ) && ~isequal( sX, sY ),
            a = false;
        elseif isa( sX{1}, 'function_handle' ),
            a = isequal( sX, sY );
        else
            for k = 1 : numel( sX ),
                [ ta, sX{k} ] = size_compat( sX{k}, sY{k} );
                a = a && ta;
            end
        end
    otherwise,
        a = isequal( sX, sY );
end
if ~a,
    sX = [];
end

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.

