function prob = ddaecoll_sol2seg(prob, oid, varargin)
% Copyright (C) Zaid Ahsan, Mingwu Li

tbid = coco_get_id(oid, 'ddaecoll'); % Create toolbox instance identifier
str  = coco_stream(varargin{:}); % Convert varargin to stream of tokens for argument parsing
run  = str.get;
if ischar(str.peek)
  soid = str.get;
else
  soid = oid;
end
lab = str.get;

[sol data] = ddaecoll_read_solution(soid, run, lab);  % Extract solution and toolbox data from disk
data       = ddaecoll_get_settings(prob, tbid, data); % Get toolbox settings
data       = ddaecoll_init_data(data, sol.x, sol.y, sol.p);  % Build toolbox data
sol        = ddaecoll_init_sol(data, sol.t, sol.x, sol.y, sol.p);  % Build initial solution guess
prob       = ddaecoll_construct_seg(prob, tbid, data, sol); % Append continuation problem

end
