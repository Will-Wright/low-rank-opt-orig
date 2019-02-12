function y = norm(v,type)

switch lower(type)
    case 'fro'
        y = sqrt(tfocs_normsq(v));
    otherwise
        error('this norm is not supported for packSVD objects');
end
