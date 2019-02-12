function [ v, g ] = AugNegLogLikelihood(U,Af,b,type,lambda,W)

%  [ v, g ] = AugNegLogLikelihood(U,Af,y,type,lambda,W)
% Computes value and gradient (at X) of augmented negative log likelihood
% 
% l(X) + lambda <W,X>
%
% where 
%
% l(X) = 0.5 ||A(X) - b||^2 (Gaussian) 
%
% l(X) = sum_k <A_k, X> - b_k log <A_k, X> (Poisson)


if strcmpi(type, 'gaussian') 
    r = Af(U) - b; 
    v = 0.5 * tfocs_normsq( r(:) ) + lambda * innprodW(U,W);
    g = r;
elseif strcmpi(type, 'poisson')
    Afx = Af(U); 
    v = sum(Afx(:)) - b(:)'* log(Afx(:)) + lambda * innprodW(U,W);
    g = 1 - b./Afx;
else error('Likelihood type is unrecognized') 
end
   
