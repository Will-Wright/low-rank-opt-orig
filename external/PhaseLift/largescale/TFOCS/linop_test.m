function varargout = linop_test( op, cmode, maxits )

%LINOP_TEST Performs an adjoint test on a linear operator.
%    LINOP_TEST( OP ) attempts to verify that a linear operator OP obeys
%    the inner product test: <A*x,y> = <x,A'*y> for all x, y. OP must be a
%    TFOCS linear operator with hard-coded size information; that is,
%    OP([],0) must return valid size info.
%
%    When called with a single argument, LINOP_TEST creates real test
%    vectors for X and Y. To test complex operators, use the two-argument
%    version LINOP_TEST( OP, cmode ), where:
%        cmode = 'R2R': real input, real output
%        cmode = 'R2C': real input, complex output
%        cmode = 'C2R': imag input, imag output
%        cmode = 'C2C': complex input, complex output
%
%    LINOP_TEST( OP, CMODE, MAXITS ) performs MAXITS iterations of the
%    test loop. MAXITS=25 is the default.
%
%   NORM = LINOP_TEST(...) returns an estimate of the norm of
%       the linear operator.

error(nargchk(1,3,nargin));
if isnumeric( op ),
    op = linop_matrix( op, 'C2C' );
end
if nargin < 2 || isempty( cmode ),
    x_real = true;
    y_real = true;
else
    switch upper( cmode ),
        case 'R2R', x_real = true; y_real = true;
        case 'R2C', x_real = true; y_real = false;
        case 'C2R', x_real = false; y_real = true;
        case 'C2C', x_real = false; y_real = false;
        otherwise, error( 'Invalid cmode: %s', cmode );
    end
end
if nargin < 3 || isempty(maxits),
    maxits = 25;
end
sz = op([],0);
if ~iscell(sz)
    sz = { [sz(2),1], [sz(1),1] };
end
nf = 0;
na = 0; 
errs = zeros(1,maxits+1);
nxe = 0; nye = 0;
for k = 1 : maxits,
    
    %
    % The adjoint test
    %
    
    if x_real,
        x = randn(sz{1});
    else
        x = randn(sz{1})+1j*randn(sz{1});
    end
    
    if y_real,
        y = randn(sz{2});
    else
        y = randn(sz{2})+1j*randn(sz{2});
    end
    
    nx = norm(x);
    Ax = op(x,1);
    nf = max( nf, norm(Ax)/nx );
    Ax_y = tfocs_dot( Ax, y ); 
    
    ny = norm(y);
    Ay = op(y,2);
    na = max( na, norm(Ay) / ny );
    Ay_x = tfocs_dot( x, Ay ); 
    
    errs(k) = abs(Ax_y-Ay_x)/(nx*ny);
    
    %
    % The norm iteration
    %
    
    if nxe == 0,
        if x_real,
            xx = randn(sz{1});
        else
            xx = randn(sz{1}) + 1j*randn(sz{1});
        end
        nxe = norm(xx);
    end
    yy = op(xx/nxe,1);
    nye = max(realmin,norm(yy));
    xx = op(yy/nye,2);
    nxe = norm(xx);
    
end

%
% Use the estimated singular vectors for a final adjoint est
%

if nxe > 0,
    Ax_y = tfocs_dot( op(xx,1), yy );
    Ay_x = tfocs_dot( op(yy,2), xx );
    errs(end) = abs(Ax_y-Ay_x) / (nxe*nye);
end

%
% Display the output
% 

nmax = max(nye,nxe);
norm_err = abs(nye-nxe) / nmax;
peak_err = max(errs) / nmax;
mean_err = mean(errs) / nmax;
rc = { 'complex', 'real' };
fprintf( 'TFOCS linear operator test:\n' );
fprintf( '   Input size:  [' ); fprintf( ' %d', sz{1} ); fprintf( ' ], %s\n', rc{x_real+1} );
fprintf( '   Output size: [' ); fprintf( ' %d', sz{2} ); fprintf( ' ], %s\n', rc{y_real+1} );
fprintf( 'After %d iterations:\n', maxits  );
fprintf( '    Norm estimates (forward/adjoint/error): %g/%g/%g\n', nye, nxe, norm_err );
fprintf( '       Gains: forward %g, adjoint %g\n', nf, na );
fprintf( '    Inner product error:\n' );
fprintf( '       Mean (absolute/relative): %g/%g\n', mean(errs), mean_err );
fprintf( '       Peak (absolute/relative): %g/%g\n', max(errs), peak_err );

if nargout > 0
    varargout{1} = mean([nye,nxe]);
end

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
