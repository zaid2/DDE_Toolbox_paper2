function prob = adjt_msddaebvp_sol2segs(prob, oid, varargin)
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

if ~isempty(str.peek)
  pt = str.get;
else
  pt = [];
end

data = coco_read_solution(ttbid, run, lab);
for i=1:data.nsegs
  toid = coco_get_id(ttbid, sprintf('seg%d', i));  % Target 'coll' object instance identifier
  soid = coco_get_id(stbid, sprintf('seg%d', i));  % Source 'coll' object instance identifier
  prob = adjt_ddaecoll_sol2seg(prob, toid, run, soid, lab, pt); % Construct 'coll' instance
end



[sol, ~] = msddaebvp_read_adjoint(oid,run,lab,pt);
data = msddaebvp_adjt_init_data(prob, ttbid, data);  % Build toolbox data
prob = msddaebvp_close_segs_adjt(prob, ttbid, data, sol); % Append boundary conditions

end
