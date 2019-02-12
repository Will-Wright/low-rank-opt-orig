% TFOCS_BACKTRACK
% Backtracking helper script.

% Quick exit for no backtracking
if beta >= 1, break; end % SRB changing == to >=

% Quick exit if no progress made
%xy = x - y;
%xy_sq = tfocs_normsq( xy );

nxx = norm(x'*x,'fro')^2;
nyy = norm(y'*y,'fro')^2;

xy_sq = max(0, nxx + nyy - 2*norm(y'*x,'fro')^2);

%if xy_sq == 0, localL = Inf; break; end
if xy_sq/sqrt(nxx*nyy) < 100*eps, localL = Inf; break; end
Afx = Af(x);
Afy = Af(y);

% Compute Lipschitz estimate
if backtrack_simple,
    if isinf( f_x ),
        f_x = apply_smooth( A_x );
    end
    Wxy = innprodW(x,W) - innprodW(y,W);
    q_x = f_y + tfocs_dot( Afx - Afy, g_y ) + lambda*Wxy + 0.5 * L * xy_sq;
    localL = L + 2 * max( f_x - q_x, 0 ) / xy_sq;
    backtrack_simple = abs( f_y - f_x ) >= backtrack_tol * max( abs( f_x ), abs( f_y ) );
%localL
else
    if isempty( g_Ax ),
        [ f_x, g_Ax ] = apply_smooth( A_x );
    end
    localL = 2 * tfocs_dot( Afx - Afy, g_Ax - g_Ay ) / xy_sq;
%disp(['non-simple ' num2str(localL)]);
end


% Exit if Lipschitz criterion satisfied, or if we hit Lexact
backtrack_steps = backtrack_steps + 1;
if localL <= L || L >= Lexact, break; end
L = min( Lexact, max( localL, L / beta ) );

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
