function prob = adjt_msddaebvp_isol2segs(prob, oid, varargin)
% Copyright (C) Zaid Ahsan

tbid = coco_get_id(oid, 'msbvp'); % Create toolbox instance identifier
data = coco_get_func_data(prob,tbid,'data');


str  = coco_stream(varargin{:});  % Convert varargin to stream of tokens for argument parsing

for i=1:data.nsegs
 segoid = coco_get_id(tbid, sprintf('seg%d', i)); % Create unique object instance identifier
 prob   = adjt_ddaecoll_isol2seg(prob, segoid, str); % Use 'coll' constructor to parse one instance of COLL
end



% msddebvp_arg_check(prob, tbid, data);         % Validate input

[sol,~] = msddaebvp_read_adjoint('','',data);
data = msddaebvp_adjt_init_data(prob, tbid, data);  % Build toolbox data
prob = msddaebvp_close_segs_adjt(prob, tbid, data, sol); % Append boundary conditions

end

function flag = is_empty_or_func(x)
%IS_EMPTY_OR_FUNC   Check if input is empty or contains a function handle.
flag = isempty(x) || isa(x, 'function_handle');
end
