function prob = adjt_msddaebvp_sol2segs(prob, oid, varargin)
%MSBVP_SOL2SEGS   Append 'msbvp' instance constructed from saved data.
%
% Support restarting continuation from a previously obtained solution,
% stored to disk.
%
% PROB     = MSBVP_SOL2SEGS(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB }
%
% PROB - Continuation problem structure.
% OID  - Target object instance identifier (string).
% RUN  - Run identifier (string).
% SOID - Source object instance identifier (string).
% LAB  - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: msbvp_sol2segs.m 2839 2015-03-05 17:09:01Z fschild $

ttbid = coco_get_id(oid, 'msbvp'); % Create toolbox instance identifier
str = coco_stream(varargin{:});    % Convert varargin to stream of tokens for argument parsing
run = str.get;
if ischar(str.peek)
  stbid = coco_get_id(str.get, 'msbvp');
else
  stbid = ttbid;
end
lab = str.get;

if ~isempty(str.peek)
  pt = str.get;
else
  pt = [];
end

data = coco_read_solution(ttbid, run, lab);
for i=1:data.nsegs
  toid = coco_get_id(ttbid, sprintf('seg%d', i));  % Target 'coll' object instance identifier
  soid = coco_get_id(stbid, sprintf('seg%d', i));  % Source 'coll' object instance identifier
  prob = adjt_ddaecoll_sol2seg(prob, toid, run, soid, lab, pt); % Construct 'coll' instance
end



[sol, ~] = msddaebvp_read_adjoint(oid,run,lab,pt);
data = msddaebvp_adjt_init_data(prob, ttbid, data);  % Build toolbox data
prob = msddaebvp_close_segs_adjt(prob, ttbid, data, sol); % Append boundary conditions

end
