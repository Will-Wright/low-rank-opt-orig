%SOLVER_ITERATE    TFOCS helper script
%	Performs the iterate processing common to all of the first-order solvers
%
%   Major inputs: x, x_old, xy , A_x, A_y, g_Ax, g_Ay 
%       (does not need g_x or g_y really)


% Test for positive stopping criteria
n_iter = n_iter + 1;
xx = x'*x;
norm_x = sqrt( tfocs_normsq( xx ) );
if size(x_old,2) ~= size(x,2); 
   x_old = x_old(:,1:size(x,2));
end
norm_dx = sqrt( tfocs_normsq( xx - x_old'*x_old ) );
%norm_x = sqrt( tfocs_normsq( x*x' ) );
%norm_dx = sqrt( tfocs_normsq( x*x' - x_old*x_old' ) );
% We also handle legacy stopping criteria used in the paper:
if stopCrit == 2 && beta >= 1
    % xy_sq not already computed
    nxx = norm(x'*x,'fro')^2;
    nyy = norm(y'*y,'fro')^2;
    xy_sq = max(nxx + nyy - 2*norm(y'*x,'fro')^2 , 0);

    %xy      = x - y;
    %xy_sq   = tfocs_normsq( xy );
end
if isnan( f_y ),
	status = 'NaN found -- aborting';
elseif (stopCrit == 1) && norm_dx  == 0
    if n_iter > 1
        status = 'Step size tolerance reached (||dx||=0)';
    end
elseif (stopCrit == 1) && norm_dx < tol * max( norm_x, 1 ),
    status = 'Step size tolerance reached';
elseif (stopCrit == 2) && L*sqrt(xy_sq) < tol,
    status = 'Step size tolerance reached';
elseif n_iter == maxIts,
	status = 'Iteration limit reached';
elseif countOps && any( tfocs_count___ >= maxCounts ),
	status = 'Function/operator count limit reached';
elseif backtrack_steps > 0 && xy_sq == 0,
    status = 'Unexpectedly small stepsize';
end

%
% For stopping criteria or algorithm control, we assume that the user
% needs the objective function value, but does not wish to do any more
% computation than necessary. So we will use the function value for y
% instead of x if that is cheaper to obtain. So here we determine what
% the additional computational costs will be, and choose the path that
% minimizes them, favoring x in the case of a tie.
%

if isempty(status) && ( ~isempty(stopFcn) || restart < 0 ),
    comp_x = [ isinf(f_x), saddle*~isempty(stopFcn)*isempty(g_Ax), isinf(C_x) ];
    comp_y = [ isinf(f_y), saddle*~isempty(stopFcn)*isempty(g_Ay), isinf(C_y) ];
    if sum(comp_x) <= sum(comp_y),
        if comp_x(2), [f_x,g_Ax] = apply_smooth(A_x);
        elseif comp_x(1), f_x = apply_smooth(A_x); end
        %if comp_x(3), C_x = proj_tsvd(x,RNK); end
        cur_pri = x; cur_dual = g_Ax;
        f_v = maxmin*(f_x+C_x);
    else
        if comp_y(2), [f_y,g_Ay] = apply_smooth(A_y);
        elseif comp_y(1), f_y = apply_smooth(A_y); end
        %if comp_y(3), C_y = proj_tsvd(y,RNK); end
        cur_pri = y; cur_dual = g_Ay;
        f_v = maxmin*(f_y+C_y);
    end
    for err_j = 1 : numel(stopFcn),
        if saddle,
            stop = stopFcn{err_j}(f_v,cur_pri, get(cur_dual, saddle_ndxs) );
        else
            stop = stopFcn{err_j}(f_v,cur_pri);
        end
        if stop
            status = sprintf('Reached user''s supplied stopping criteria no. %d',err_j);
        end
    end
end

%
% Data collection. Since this is used only for algorithm analysis and not
% for algorithm control, we are free to make additional computations here
% without adding them to our total algorithm cost. We prefer to track the
% x sequence in this case, since it is what we will ultimately choose as
% the solution at the end of the algorithm.
%

will_print = fid && printEvery && ( ~isempty( status ) || ~mod( n_iter, printEvery ) );
if saveHist || will_print,
    f_x_save = f_x;
    g_Ax_save = g_Ax;
    if ~isempty(errFcn) && saddle,
        if isempty(g_Ax),
            [ f_x, g_Ax ] = smoothF( A_x );
        end
        out.dual = get( g_Ax, saddle_ndxs );
    end
    if isinf(f_x),
        f_x = smoothF( A_x );
    end
    if isinf( C_x ),
        %C_x = proj_tsvd( x, RNK ); 
    end
    %f_w = maxmin * ( f_x + C_x );
    f_w = maxmin * f_x;
    if ~isempty(errFcn) && iscell(errFcn)
        for err_j = 1 : numel(errFcn),
            if saddle,
                errs(err_j) = errFcn{err_j}(f_w,x,out.dual);
            else
                errs(err_j) = errFcn{err_j}(f_w,x);
            end
        end
    end
    f_x = f_x_save;
    g_Ax = g_Ax_save;
end

% Register a warning if the step size suggests a Lipschitz violation
if isempty(status) && ( beta < 1 && backtrack_simple && localL > Lexact ),
    warning_lipschitz = true;
end	

% Print status
if will_print,
    if warning_lipschitz,
        warning_lipschitz = false;
        bchar = 'L';
    elseif backtrack_simple,
        bchar = ' '; 
    else
        bchar = '*'; 
    end
	fprintf( fid, '%-4d| %+12.5e  %8.2e  %8.2e%c', ...
        n_iter, f_w, norm_dx / max( norm_x, 1 ), 1 / L, bchar );
    if countOps,
        fprintf( fid, '|' );
        fprintf( fid, ' %5d', tfocs_count___ );
    end
    if ~isempty(errFcn),
        if countOps,
            fprintf( fid, ' ' );
        end
    	fprintf( fid, '|' );
        fprintf( fid, ' %8.2e', errs );
    end
	fprintf( fid, '\n');
end

% Save history, extending arrays if needed
if saveHist,
    if length(out.f) < n_iter && isempty(status),
        csize = min(maxIts,length(out.f)+1000);
        out.f(end+1:csize,1) = 0;
        out.theta(end+1:csize,1) = 0;
        out.stepsize(end+1:csize,1) = 0;
        out.normGrad(end+1:csize,1) = 0;
        if countOps,
            out.counts(end+1:csize,:) = 0;
        end
        if ~isempty(errs),
            out.err(end+1:csize,:) = 0;
        end
    end
    out.f(n_iter) = f_w;
    out.theta(n_iter) = theta;
    out.stepsize(n_iter) = 1 / L;
    out.normGrad(n_iter) = norm_dx;
    if countOps,
        out.counts(n_iter,:) = tfocs_count___;
    end
    if ~isempty(errFcn),
        out.err(n_iter,:) = errs;
    end
end

% Exit if positive stopping criteria met
if ~isempty( status ),
	break;
end

% Restart acceleration if necessary
backtrack_steps = 0;
if n_iter - restart_iter == abs(restart) || restart < 0 && maxmin * f_v > f_v_old,
disp('restart accel')

    restart_iter = n_iter;
    backtrack_simple = true;
	theta = Inf;
    y = x; A_y = A_x; f_y = f_x; g_Ay = g_Ax; g_y = g_x; C_y = C_x;
    z = x; A_z = A_x; f_z = f_x; g_Az = g_Ax; g_z = g_x; C_z = C_x;
    f_v_old = Inf;
    continue;
elseif restart < 0,
    f_v_old = f_v;
end

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.


