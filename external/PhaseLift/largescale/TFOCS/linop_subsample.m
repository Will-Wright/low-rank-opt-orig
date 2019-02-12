function op = linop_subsample( sz, omega )

%OP = LINOP_SUBSAMPLE( SZ, OMEGA )

% Designed to be used with a fft or dct
% Do we need a separate oversampling version?
% Also, make linop_DCT and linop_FFT ?
%   The reason we might want this is that for linop_FFT,
%   people will probably use idct as the transpose -- this is not 
%   correct, due to scaling

% Help documentation: TBD
%    Constructs a TFOCS-compatible linear operator from separate function
%    handles that compute the forward and adjoint operations. The first
%    argument, SZ, gives the size of the linear operator; and the forward
%    and adjoint handles are AF and AT, respectively.
% 
%    If the inputs and outputs are simple vectors, then SZ can take the
%    standard Matlab form [N,M], where N is the input length and M is the
%    output length. If the input or output is a matrix or array, then SZ
%    should take the form { S_in, S_out }, where S_in is the size of the
%    input and S_out is the size of the output.
%

error(nargchk(1,2,nargin));
if numel( sz ) ~= 2,
    error( 'Size must have two elements.' );
elseif isnumeric( sz ),
    sz = { [sz(2),1], [sz(1),1] };
elseif ~isa( sz, 'cell' ),
    error( 'Invalid operator size specification.' );
end

dom = sz{1};
n1 = dom(1); n2 = dom(2);
ran = sz{2};

% There are several possibilities.  Let x be a point in the domain
%{
    x is a scalar.  Then omega should be a vector.
                    We return x(omega), resp.

    x is a matrix.
        omega is a vector
            This is ambiguous: does the user want x(omega,:) or x(omega)?
        omega is a matrix with 2 columns, then assume it is [I,J]
            We convert it to linear indices.
        omega is a general matrix, more than 2 columns
            Not sure what the user means; report an error.
        omega is a sparse matrix
            We find it's entries, and use those

%}

if n2 == 1
    % x is a vector, not a matrix.  Simple.
    op = @(x,mode) linop_subsample_vector( sz, omega, x, mode );
else
    % trickier case.
    if issparse(omega)
        indx = find(omega);
        op = @(x,mode) linop_subsample_vector( sz, indx, x, mode );
    elseif isvector(omega)
        op = @(x,mode) linop_subsample_vector( sz, omega, x, mode );
    elseif size(omega,2) == 2
        indx = sub2ind( sz{1}, omega(:,1), omega(:,2) );
        op = @(x,mode) linop_subsample_vector( sz, indx, x, mode );
    else
        error('omega is not an acceptable size');
    end
        
    
end


function y = linop_subsample_vector(sz, omega, x, mode )
switch mode,
    case 0, y = sz;
    case 1, y = x(omega);
    case 2, 
        y = zeros( sz{1} );
        y(omega) = x;
end



% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
