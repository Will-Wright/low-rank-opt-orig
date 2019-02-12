function re = isreal(v)

re = true;
if ~isempty(v.X)
    re = re & isreal(v.X);
end
if ~isempty(v.s)
    re = re & isreal(v.s) & isreal(v.U) & isreal(v.V);
end