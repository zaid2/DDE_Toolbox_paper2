function prob = adjt_ddaecoll_sol2seg(prob, oid, varargin)
%ADJT_COLL2COLL   Append adjoint of 'coll' instance from saved solution.
%
% PROB = ADJT_COLL2COLL(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB }
%
% Append adjoint of a 'coll' instance with object instance identifier OID
% that has been previously added to the continuation problem contained in
% PROB from the same saved solution using ODE_COLL2COLL.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string). Pass the empty string
%        '' for a simple continuation of trajectory segments.
%
% See ODE_ISOL2COLL for more details on PROB and OID.
%
% RUN  : Run identifier (string or cell-array of strings). Name of the run
%        from which to restart a new continuation run.
% SOID : Source object instance identifier (string, optional). If the
%        argument SOID is omitted, OID is used. Pass the empty string ''
%        for OID and omit SOID for a simple continuation of trajectory
%        segments. Pass non-trivial object identifiers if an instance of
%        the COLL toolbox is part of a composite continuation problem.
% LAB  : Solution label (integer). The integer label assigned by COCO to an
%        trajectory segment during the continuation run RUN.
%
% See also: ADJT_ISOL2COLL, ADJT_BP2COLL, COLL_READ_ADJOINT,
% COLL_ADJT_INIT_DATA, COLL_CONSTRUCT_ADJT

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: ode_coll2coll.m 2898 2015-10-07 21:17:13Z hdankowicz $

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
