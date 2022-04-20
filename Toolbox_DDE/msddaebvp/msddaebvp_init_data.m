function data = msddaebvp_init_data(prob, tbid, data)
% Copyright (C) Zaid Ahsan

xnum = 0;
for i=1:data.nsegs
  stbid = coco_get_id(tbid, sprintf('seg%d.ddaecoll', i)); % Construct 'coll' toolbox instance identifier
  fdata = coco_get_func_data(prob, stbid, 'data');     % Extract 'coll' toolbox data
  xnum  = xnum+numel(fdata.x0_idx);                    % Track total dimension
end

data.T0_idx = (1:data.nsegs)';                     % Index array for initial times
data.T_idx  = data.nsegs+(1:data.nsegs)';          % Index array for interval lengths
data.x0_idx = 2*data.nsegs+(1:xnum)';              % Index array for trajectory end points at t=0
data.x1_idx = 2*data.nsegs+xnum +(1:xnum)';        % Index array for trajectory end points at t=1
data.p_idx  = 2*data.nsegs+2*xnum+(1:fdata.pdim)'; % Index array for problem parameters

end
