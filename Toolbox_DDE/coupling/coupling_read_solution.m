function [sol, data] = coupling_read_solution(tbid, run, lab)

% tbid         = coco_get_id(oid, 'coupling');
[data,chart] = coco_read_solution(tbid, run, lab);


sol.Delta = chart.x(data.Delta_idx);
sol.gamma = chart.x(data.gamma_idx);
if data.pnew_dim~=0
sol.p0    = chart.x(data.pnew_idx);
else
    sol.p0 = [];
end

end