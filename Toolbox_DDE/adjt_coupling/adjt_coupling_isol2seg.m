function prob = adjt_coupling_isol2seg(prob, oid, W, varargin)
%ADJT_ISOL2COLL   Append adjoint of 'coll' instance from initial guess.
%
% PROB = ADJT_ISOL2COLL(PROB, OID, VARARGIN)
% VARARGIN = { [OPTS] }
% OPTS = { '-coll-end' | '-end-coll' }
%
% Append adjoint of a 'coll' instance with object instance identifier OID
% that has been previously added to the continuation problem contained in
% PROB using ODE_ISOL2COLL. The preceding call to ODE_ISOL2COLL must
% include explicit Jacobians, while functions evaluating second derivatives
% are optional. 
%
% On input:
%
% PROB   : Continuation problem structure.
%
% OID    : Object instance identifier (string). The corresponding toolbox
%          instance identifier is coco_get_id(OID, 'coll'). Pass the empty
%          string '' for a simple continuation of trajectory segments. Pass
%          a non-trivial object identifier if an instance of the COLL
%          toolbox is part of a composite continuation problem.
%
% OPTS   : '-coll-end' and '-end-coll' (optional, multiple options may be
%          given). Either '-coll-end' or '-end-coll' marks the end of input
%          to ADJT_ISOL2COLL.
%
% See also: ODE_ISOL2COLL, COLL_READ_ADJOINT, COLL_ADJT_INIT_DATA,
% COLL_CONSTRUCT_ADJT

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: ode_isol2coll.m 2898 2015-10-07 21:17:13Z hdankowicz $

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








