function prob = adjt_msddaebvp_isol2segs(prob, oid, varargin)
%MSBVP_ISOL2SEGS   Append 'msbvp' instance constructed from initial data.
%
% Construct multiple instances of 'coll' and append boundary conditions.
%
% PROB     = MSBVP_ISOL2SEGS(PROB, OID, VARARGIN)
% VARARGIN = { COLL... ( PNAMES | 'END-COLL' ) BCND }
% BCND     = @BC [@DBCDX] [BC_DATA [@BC_UPDATE]]
%
% PROB       - Continuation problem structure.
% OID        - Object instance identifier (string).
% COLL       - Argument to coll_isol2seg without PNAMES.
% PNAMES     - Optional string label or cell array of string labels for
%              continuation parameters tracking problem parameters.
% @BC        - Boundary conditions function handle.
% @DBCDX     - Optional function handle to Jacobian w.r.t. T, x0, x1, p.
% BC_DATA    - Optional boundary condition function data (struct).
% @BC_UPDATE - Optional function handle to function data updater.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: msbvp_isol2segs.m 2839 2015-03-05 17:09:01Z fschild $

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
