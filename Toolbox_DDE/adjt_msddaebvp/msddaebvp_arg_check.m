function msddebvp_arg_check(prob, tbid, data)
% Copyright (C) Zaid Ahsan

assert(data.nsegs~=0, '%s: insufficient number of segments', tbid);
pnum = [];
for i=1:data.nsegs
  fid   = coco_get_id(tbid,sprintf('seg%d.ddecoll', i));
  fdata = coco_get_func_data(prob, fid, 'data');
  assert(isempty(fdata.pnames), ...
    '%s: parameter labels must not be passed to coll', tbid);
  assert(isempty(pnum) || pnum==numel(fdata.p_idx), '%s: %s', ...
    tbid, 'number of parameters must be equal for all segments');
  pnum = numel(fdata.p_idx);
end
assert(iscellstr(data.pnames) || isempty(data.pnames), ...
  '%s: incorrect format for parameter labels', tbid);
assert(pnum==numel(data.pnames) || isempty(data.pnames), ...
  '%s: incompatible number of parameter labels', ...
  tbid);
assert(isa(data.fbchan, 'function_handle'), ...
  '%s: input for ''fbc'' is not a function handle', tbid);

end
