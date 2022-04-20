function prob = adjt_ddaecoll_sol2seg(prob, oid, varargin)
% Copyright (C) Zaid Ahsan, Mingwu Li

grammar   = 'RUN [SOID] LAB [POINT] [OPTS]';
args_spec = {
   'RUN', 'cell', '{str}',   'run',  {}, 'read', {}
  'SOID',     '',   'str',  'soid', oid, 'read', {}
   'LAB',     '',   'num',   'lab',  [], 'read', {}
 'POINT',     '',   'str', 'point',  {}, 'read', {}
  };
[args, opts] = coco_parse(grammar, args_spec, [], varargin{:});


tbid = coco_get_id(oid, 'ddaecoll');
data = coco_get_func_data(prob, tbid, 'data');


sol  = ddaecoll_read_adjoint(args.soid, args.run, args.lab, args.point);
data = ddaecoll_adjt_init_data(prob, data);
prob = ddaecoll_construct_adjt(prob, tbid, data, sol);

end
