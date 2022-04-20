function prob = msddaebvp_sol2segs(prob, oid, varargin)
% Copyright (C) Zaid Ahsan

ttbid = coco_get_id(oid, 'msbvp'); % Create toolbox instance identifier
str = coco_stream(varargin{:});    % Convert varargin to stream of tokens for argument parsing
run = str.get;
if ischar(str.peek)
  stbid = coco_get_id(str.get, 'msbvp');
else
  stbid = ttbid;
end
lab = str.get;

data = coco_read_solution(ttbid, run, lab);
for i=1:data.nsegs
  toid = coco_get_id(ttbid, sprintf('seg%d', i));  % Target 'coll' object instance identifier
  soid = coco_get_id(stbid, sprintf('seg%d', i));  % Source 'coll' object instance identifier
  prob = ddaecoll_sol2seg(prob, toid, run, soid, lab); % Construct 'coll' instance
end
data = msddaebvp_init_data(prob, ttbid, data);  % Build toolbox data
prob = msddaebvp_close_segs(prob, ttbid, data); % Append boundary conditions

end
