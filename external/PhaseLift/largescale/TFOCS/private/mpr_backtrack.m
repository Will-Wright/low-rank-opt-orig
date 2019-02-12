% TFOCS_BACKTRACK
% Backtracking helper script.

% Quick exit for no backtracking
if beta >= 1, break; end % SRB changing == to >=

% Quick exit if no progress made
%size(x),size(y)
xy = x - y;
%xy_sq = tfocs_normsq( xy );
prod_yx = y'*x;

xy_sq = norm(x'*x,'fro')^2 + norm(y'*y,'fro')^2 - 2*norm(prod_yx,'fro')^2;

if xy_sq == 0, localL = Inf; break; end

% Compute Lipschitz estimate
if backtrack_simple,
    if isinf( f_x ),
        f_x = apply_smooth( A_x );
    end
    %q_x = f_y + tfocs_dot( Af(x)-Af(y), g_y ) + 0.5 * L * xy_sq;
    if epsilon == Inf
       Wxy = tfocs_dot(x,x) - tfocs_dot(y,y);
    else
       Wxy = tfocs_dot(x,x)-tfocs_dot(x'*W)+tfocs_dot(y'*W)-tfocs_dot(y,y);
       Wxy = 1/epsilon*Wxy;
    end
    q_x = f_y + tfocs_dot( Af(x)-Af(y), g_y ) + lambda*Wxy + 0.5 * L * xy_sq;
    localL = L + 2 * max( f_x - q_x, 0 ) / xy_sq;
    backtrack_simple = abs( f_y - f_x ) >= backtrack_tol * max( abs( f_x ), abs( f_y ) );
backtrack_simple = 1;
else
    if isempty( g_Ax ),
        [ f_x, g_Ax ] = apply_smooth( A_x );
    end
    localL = 2 * tfocs_dot( A_x - A_y, g_Ax - g_Ay ) / xy_sq;
end

% Exit if Lipschitz criterion satisfied, or if we hit Lexact
backtrack_steps = backtrack_steps + 1;
if localL <= L || L >= Lexact, break; end
L = min( Lexact, max( localL, L / beta ) );

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
