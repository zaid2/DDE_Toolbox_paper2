function prob = adjt_coupling_sol2seg(prob, oid, W, varargin)

% Copyright (C) Zaid Ahsan

grammar   = 'RUN [SOID] LAB [POINT] [OPTS]';
args_spec = {
   'RUN', 'cell', '{str}',   'run',  {}, 'read', {}
  'SOID',     '',   'str',  'soid', oid, 'read', {}
   'LAB',     '',   'num',   'lab',  [], 'read', {}
 'POINT',     '',   'str', 'point',  {}, 'read', {}
  };
[args, opts] = coco_parse(grammar, args_spec, [], varargin{:});

% if opts.switch
%   prob = adjt_BP2coll(prob, oid, args.run, args.soid, args.lab);
%   return
% end


% fid  = sprintf('gseg.%s',gidx);
% soid = coco_get_id(args.soid, 'coup_dde');
% soid = coco_get_id(soid, fid);
% data = coco_get_func_data(prob, soid, 'data');
% 
% 
% tbid      = coco_get_id(oid, 'coup_dde');
% data.gidx = gidx;

% data.axidx = [];
% for i=1:data.ncoup
% fname = sprintf('seg%d',i);          
% fbid              = coco_get_id(fname, 'ddaecoll'); % Create toolbox instance identifier
% [fdata, axidx]    = coco_get_adjt_data(prob, fbid, 'data', 'axidx');
% data.axidx = [data.axidx; axidx];
% end
% data.ddaecoll_opt = fdata.ddaecoll_opt;


if isempty(W)
tbid = coco_get_id(oid, 'coupling'); % Create toolbox instance identifier
else
fname = sprintf('coupling%d',W);
tbid = coco_get_id(oid,fname);
end
data = coco_get_func_data(prob,tbid,'data');

fbid     = coco_get_id(oid, 'ddaecoll');
fdata    = coco_get_adjt_data(prob, fbid, 'data');
data.ddaecoll_opt = fdata.ddaecoll_opt;

% sol  = coupling_read_adjoint(args.soid, gidx, args.run, args.lab, args.point);
[sol,~] = coupling_read_adjoint(oid, W, args.run, args.lab, args.point);
data = adjt_coupling_init_data(prob, data);
prob = coupling_construct_adjt(prob, tbid, data, sol);







end















