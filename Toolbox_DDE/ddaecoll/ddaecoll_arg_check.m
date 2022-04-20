function ddaecoll_arg_check(tbid, data, t0, x0, y0, p0)
% Copyright (C) Zaid Ahsan, Mingwu Li

assert(isa(data.fhan, 'function_handle'), ...
  '%s: input for ''f'' is not a function handle', tbid);
assert(isnumeric(t0), '%s: input for ''t0'' is not numeric', tbid);
assert(isnumeric(x0), '%s: input for ''x0'' is not numeric', tbid);
assert(isnumeric(y0), '%s: input for ''y0'' is not numeric', tbid);
assert(isnumeric(p0), '%s: input for ''p0'' is not numeric', tbid);
assert(ndims(t0)==2 && min(size(t0))==1, ...
  '%s: input for ''t0'' is not a vector', tbid);
assert(size(t0,1)>=size(t0,2),...
  '%s: input for ''t0'' is not a column vector', tbid); 
assert(ndims(x0)==2, ...
  '%s: input for ''x0'' is not an array of vectors', tbid);
assert(size(x0, 1)==numel(t0), ...
  '%s: dimensions of ''t0'' and ''x0'' do not match', tbid);
assert(ndims(y0)==2, ...
  '%s: input for ''y0'' is not an array of vectors', tbid);
if ~isempty(y0)
assert(size(y0, 1)==numel(t0), ...
  '%s: dimensions of ''t0'' and ''y0'' do not match', tbid);
end
assert(numel(p0)==numel(data.pnames) || isempty(data.pnames), ...
  '%s: incompatible number of elements for ''p0'' and ''pnames''', ...
  tbid);



end
