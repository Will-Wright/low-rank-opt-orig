function ans = size_ambig( sz )

switch class( sz ),
    case 'double',
        ans = isempty( sz );
    case 'cell',
        if isa( sz{1}, 'function_handle' ),
            ans = false;
        elseif isempty( sz ),
            ans = true;
        else
            ans = false;
            for k = 1 : numel( sz ),
                if size_ambig( sz{k} ),
                    ans = true;
                    break;
                end
            end
        end
    otherwise,
        ans = true;
end

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
