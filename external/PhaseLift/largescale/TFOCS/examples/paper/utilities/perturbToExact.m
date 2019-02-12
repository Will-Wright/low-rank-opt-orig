function [xNew,bNew,Dnew,lambda,z] = perturbToExact(lambda,x,Af,At,D,delta0,b,mu,xPlug)
%  [x,b,D] = perturbToExact(lambda,x,Af,At,D,delta0,b,mu,xPlug)
%      takes a near Dantzig solution and modifies
%      x, b and D just a bit to get an exact solution
%
%  [x,b,D,lambda,z] = ...
%       also returns the hard-thresholded dual solution lambda,
%       and the variable z := b - A(x)
%
%   TFOCS ver

% load temp



z       = b - Af(x);
Atz     = At(z);
DAtz    = D.*Atz;
fprintf('  Initial x from solver is infeasible by %.2e\n',norm(DAtz,Inf) - delta0 );

%    Cleanup the dual solution by hard thresholding
% %%
% figure(2); clf;
% [lambdaAbs,indx] = sort( abs(lambda),'descend' );
% semilogy( lambdaAbs,'o:' )
% hold all
% semilogy( delta0 - abs(DAtz(indx))   ); % there are some negatives too, from the infeasible parts
% semilogy(-(delta0 - abs(DAtz(indx)))   ); % there are some negatives too, from the infeasible parts
% %%
S       = find( abs(lambda) > 1e-5*norm(lambda,Inf) );
% S       = find(   (abs(lambda) > delta0 - abs(DAtz))   &  lambda  );
N       = length(lambda);
Sc      = setdiff(1:N,S);    % the complement of S
lambda(Sc) = 0;

if delta0 == 0
    % special case, and very easy to deal with.  The dual variable
    %   is now unconstrained, and the D parameter is irrelevant
    Dnew        = D;  % doesn't matter
    xNew        = shrink( xPlug - At(Af(Dnew.*lambda))/mu, 1/mu );
    bNew        = Af(xNew);     % so xNew is feasible by construction
    z           = zeros(size(bNew));
    return;
end
    

% on S, we want DAtz = -delta0*sign(lambda)
Dnew    = D;
Dnew(S) = -delta0*sign(lambda(S))./Atz(S);
if any( sign(Dnew(S)) ~= sign(D(S)) )
    disp('Warning: D_S has changed signs');
    % but is this bad?
end

% on S^c, don't need to modify DAtz unless it's infeasible
if norm( DAtz(Sc), Inf ) > delta0
    disp('Warning: D_Sc needs to be shrunk');
    infeas_set      = intersect( Sc, find(abs(DAtz)>delta0) );
    Dnew(infeas_set)= .99*delta0./D(infeas_set);
end

xNew        = shrink( xPlug - At(Af(Dnew.*lambda))/mu, 1/mu );
bNew        = z + Af( xNew );

%-- Note that xNew is feasible by construction, since bNew - Af(xNew) = z
% infeasibilty = norm( Dnew.*At( bNew - Af(xNew)), Inf ) - delta0;
% should be machine precision

% Now, can we verify if we have a dual solution?
%{
Solve:  sgn( x_T ) = A_T^* A_S D_S lambda_S

    May not be solvable, since |T| may be greater than |S|

    Even if it is solvable, need to make sure subgradient fits for x_T^c

I think I need to solve the auxiliary equation first, so hold off
%}
% T       = find(xNew );