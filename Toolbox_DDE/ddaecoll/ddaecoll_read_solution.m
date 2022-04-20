function [sol data] = ddaecoll_read_solution(oid, run, lab)
% Copyright (C) Zaid Ahsan, Mingwu Li

tbid         = coco_get_id(oid, 'ddaecoll');
[data chart] = coco_read_solution(tbid, run, lab);

sol.t  = chart.x(data.T0_idx)+data.tbp(data.tbp_idx)*chart.x(data.T_idx);
xbp    = reshape(chart.x(data.xbp_idx), data.xbp_shp)';
sol.x  = xbp(data.tbp_idx,:);
if ~isempty(data.ybp_idx)
ybp    = reshape(chart.x(data.ybp_idx), data.ybp_shp)';
sol.y  = ybp(data.tbp_idx,:);
else
    sol.y = [];
end
sol.p  = chart.x(data.p_idx);


end
