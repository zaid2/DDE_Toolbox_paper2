function prob = adjt_coupling_isol2seg(prob, oid, W, varargin)

% Copyright (C) Zaid Ahsan

if isempty(W)
tbid = coco_get_id(oid, 'coupling'); % Create toolbox instance identifier
else
fname = sprintf('coupling%d',W);
tbid = coco_get_id(oid,fname);
end
data = coco_get_func_data(prob, tbid, 'data');


% data.ddaecoll_seg = fdata.ddaecoll_seg;

% data.axidx        = axidx;

% data.axidx = [];
% for i=1:data.nseg
% fname = sprintf('seg%d',i);          
% fbid              = coco_get_id(fname, 'ddaecoll'); % Create toolbox instance identifier
% [fdata, axidx]    = coco_get_adjt_data(prob, fbid, 'data', 'axidx');
% data.axidx = [data.axidx; axidx];
% end

% fname = sprintf('seg%d',i); 
fbid     = coco_get_id(oid, 'ddaecoll');
fdata    = coco_get_adjt_data(prob, fbid, 'data');
data.ddaecoll_opt = fdata.ddaecoll_opt;

data = adjt_coupling_init_data(prob, data);
[sol,~]  = coupling_read_adjoint('', W, '', data);
prob = coupling_construct_adjt(prob, tbid, data, sol);

end








