function y = ctranspose(v)
y=v;
if ~isempty(v.X)
    y.X = v.X';
end
if ~isempty(v.U)
    % The matrix is U*diag(S)*V'
    % where U, S and V could be complex.
    % We want the ' transpose, so we have
    %   Y = (V')' * diag(S) * U'
    y.U = v.V;
    y.V = v.U;
    y.s = v.s;
end
