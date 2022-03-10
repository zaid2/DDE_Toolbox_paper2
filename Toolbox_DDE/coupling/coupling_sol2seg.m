function prob = coupling_sol2seg(prob, oid, W, varargin)
%COLL_SOL2SEG   Append 'coll' instance constructed from saved data.
%
% Support restarting continuation from a previously obtained solution,
% stored to disk.
%
% PROB     = COLL_SOL2SEG(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB }
%
% PROB - Continuation problem structure.
% OID  - Target object instance identifier (string).
% RUN  - Run identifier (string).
% SOID - Source object instance identifier (string).
% LAB  - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_sol2seg.m 2839 2015-03-05 17:09:01Z fschild $

if isempty(W)
tbid = coco_get_id(oid, 'coupling'); % Create toolbox instance identifier
else
fname = sprintf('coupling%d',W);
tbid = coco_get_id(oid,fname);
end

str  = coco_stream(varargin{:}); % Convert varargin to stream of tokens for argument parsing
run  = str.get;
% if ischar(str.peek)
%   soid = str.get;
% else
%   soid = oid;
% end
lab = str.get;


[sol data] = coupling_read_solution(tbid, run, lab);  % Extract solution and toolbox data from disk

% NTST = coco_get(prob,'ddaecoll','NTST');
% NCOL = coco_get(prob,'ddaecoll','NCOL');

fid = coco_get_id(oid, 'ddaecoll');
ddaecoll_seg = coco_get_func_data(prob,fid,'data');

data.ddaecoll_seg = ddaecoll_seg;
data       = coupling_init_data(data, sol.p0);  % Build toolbox data

prob       = coupling_construct_seg(prob, tbid, data, [sol.p0(:);sol.Delta(:);sol.gamma(:)]); % Append continuation problem

end
