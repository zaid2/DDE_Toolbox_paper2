function prob = adjt_ddaecoll_isol2seg(prob, oid, varargin)
% Copyright (C) Zaid Ahsan, Mingwu Li

tbid = coco_get_id(oid, 'ddaecoll');
data = coco_get_func_data(prob, tbid, 'data');

% data = ddaecoll_get_settings(prob, tbid, data);     % Get toolbox settings
data = ddaecoll_adjt_init_data(prob, data);
sol  = ddaecoll_read_adjoint('', '', data);
prob = ddaecoll_construct_adjt(prob, tbid, data, sol);

end
