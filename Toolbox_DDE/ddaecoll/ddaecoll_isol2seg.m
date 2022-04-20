function prob = ddaecoll_isol2seg(prob, oid, varargin)
% Copyright (C) Zaid Ahsan, Mingwu Li

tbid = coco_get_id(oid, 'ddaecoll'); % Create toolbox instance identifier
str  = coco_stream(varargin{:}); % Convert varargin to stream of tokens for argument parsing
data.fhan = str.get;
data.dfdthan  = [];
data.dfdxhan  = [];
data.dfdyhan  = [];
data.dfdphan  = [];
data.dfdtdthan  = [];
data.dfdtdxhan  = [];
data.dfdtdyhan  = [];
data.dfdtdphan  = [];
data.dfdxdxhan  = [];
data.dfdxdyhan  = [];
data.dfdxdphan  = [];
data.dfdydyhan  = [];
data.dfdydphan  = [];
data.dfdpdphan  = [];

if is_empty_or_func(str.peek)
  data.dfdthan = str.get;
  if is_empty_or_func(str.peek)
    data.dfdxhan = str.get;
    if is_empty_or_func(str.peek)
        data.dfdyhan = str.get; 
        if is_empty_or_func(str.peek)
            data.dfdphan = str.get;
            if is_empty_or_func(str.peek)
                data.dfdtdthan = str.get;
                if is_empty_or_func(str.peek)
                    data.dfdtdxhan = str.get;
                    if is_empty_or_func(str.peek)
                        data.dfdtdyhan = str.get;
                        if is_empty_or_func(str.peek)
                            data.dfdtdphan = str.get;
                            if is_empty_or_func(str.peek)
                                data.dfdxdxhan = str.get;
                                if is_empty_or_func(str.peek)
                                    data.dfdxdyhan = str.get;
                                    if is_empty_or_func(str.peek)
                                        data.dfdxdphan = str.get;
                                        if is_empty_or_func(str.peek)
                                            data.dfdydyhan = str.get;
                                            if is_empty_or_func(str.peek)
                                                data.dfdydphan = str.get;
                                                if is_empty_or_func(str.peek)
                                                    data.dfdpdphan = str.get;
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
  end
end
t0 = str.get;
x0 = str.get;
y0 = str.get; % a cell array where each entry corresponds to a y
data.pnames = {};
if iscellstr(str.peek('cell'))
  data.pnames = str.get('cell');
end
p0 = str.get;

% ddaecoll_arg_check(tbid, data, t0, x0, y0, p0);     % Validate input
data = ddaecoll_get_settings(prob, tbid, data);     % Get toolbox settings
data = ddaecoll_init_data(data, x0, y0, p0);        % Build toolbox data
sol = ddaecoll_init_sol(data, t0, x0, y0, p0);     % Build initial solution guess
prob = ddaecoll_construct_seg(prob, tbid, data, sol); % Append continuation problem

end

function flag = is_empty_or_func(x)
%IS_EMPTY_OR_FUNC   Check if input is empty or contains a function handle.
flag = isempty(x) || isa(x, 'function_handle');
end
