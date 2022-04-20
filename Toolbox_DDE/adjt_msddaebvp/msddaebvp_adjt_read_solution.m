function [sol data] = msddaebvp_adjt_read_solution(oid, run, lab)
% Copyright (C) Zaid Ahsan

tbid = coco_get_id(oid, 'msbvp');
data = coco_read_solution(tbid, run, lab);

sol = cell(1, data.nsegs);
for i=1:data.nsegs
  segoid = coco_get_id(tbid, sprintf('seg%d', i)); % 'coll' object instance identifier
  sol{i} = ddaecoll_read_solution(segoid, run, lab); % Trajectory segment
end

end
