function pattern = make(name, h, w, d, varargin)
p = inputParser;
p.addParameter('normalize', false);
p.addParameter('forceregular', false);
p.parse(varargin{:});
normalize = p.Results.normalize;
forceregular = p.Results.forceregular;

switch lower(name)
   case 'binary'
      pattern = cdp.binary(h,w,d);
   case 'rademacher'
      pattern = cdp.rademacher(h,w,d);
   case 'ternary'
      pattern = cdp.ternary(h,w,d);
   case 'aternary'
      pattern = cdp.aternary(h,w,d);
   case 'quaternary'
      pattern = cdp.quaternary(h,w,d);
   case 'octanary'
      pattern = cdp.octanary(h,w,d);
   case 'aoctanary'
      pattern = cdp.aoctanary(h,w,d);
   case 'runiform'
      pattern = cdp.runiform(h,w,d);
   case 'cuniform'
      pattern = cdp.cuniform(h,w,d);
   case 'rgaussian'
      pattern = cdp.rgaussian(h,w,d);
   case 'cgaussian'
      pattern = cdp.cgaussian(h,w,d);
   otherwise
      error('[cdp.make] Pattern name ''%s'' invalid!',name);
end
if forceregular
   pattern(:,:,1) = 1;
end
if normalize
   pattern = cdp.normalize(pattern);
end

end
