function Y = full(v)

if ~isempty(v.X)
    Y = full(v.X);
else
    Y = 0;
end

if ~isempty(v.s)
    Y = Y + v.U*diag(v.s)*v.V';
end

if numel(Y) == 1
    Y = repmat(Y,v.sz(1),v.sz(2) );
end
