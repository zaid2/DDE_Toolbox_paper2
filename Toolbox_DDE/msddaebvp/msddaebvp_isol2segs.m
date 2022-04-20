function prob = msddaebvp_isol2segs(prob, oid, varargin)

% Copyright (C) Zaid Ahsan

tbid = coco_get_id(oid, 'msbvp'); % Create toolbox instance identifier
str  = coco_stream(varargin{:});  % Convert varargin to stream of tokens for argument parsing
data.nsegs = 0;
while isa(str.peek, 'function_handle')
  data.nsegs = data.nsegs+1;
  segoid = coco_get_id(tbid, sprintf('seg%d', data.nsegs)); % Create unique object instance identifier
  prob   = ddaecoll_isol2seg(prob, segoid, str); % Use 'coll' constructor to parse one instance of COLL
end
data.pnames = {};
if strcmpi(str.peek, 'end-coll') % Check for stop token
  str.skip;
else
  data.pnames = str.get('cell');
end
data.fbchan = str.get;
data.dfbcdxhan = [];
if is_empty_or_func(str.peek)
  data.dfbcdxhan = str.get;
end
data.bc_data   = struct();
data.bc_update = [];
if isstruct(str.peek)
  data.bc_data = str.get;  
  if is_empty_or_func(str.peek)
    data.bc_update = str.get;
  end
end

% msddebvp_arg_check(prob, tbid, data);         % Validate input
data = msddaebvp_init_data(prob, tbid, data);  % Build toolbox data
prob = msddaebvp_close_segs(prob, tbid, data); % Append boundary conditions

end

function flag = is_empty_or_func(x)
%IS_EMPTY_OR_FUNC   Check if input is empty or contains a function handle.
flag = isempty(x) || isa(x, 'function_handle');
end
