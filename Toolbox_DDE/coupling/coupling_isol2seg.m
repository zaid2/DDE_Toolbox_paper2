function prob = coupling_isol2seg(prob, oid, W, varargin)

if isempty(W)
tbid = coco_get_id(oid, 'coupling'); % Create toolbox instance identifier
else
fname = sprintf('coupling%d',W);
tbid = coco_get_id(oid,fname);
end
str  = coco_stream(varargin{:}); % Convert varargin to stream of tokens for argument parsing

data.xi = str.get;   %==== i^{th} segment
yiid = str.get;
data.yiid = yiid;
gseg = str.get;


Deltai = []; gammai = [];
for k=1:length(gseg)
Xj{k} = gseg{k}{1};     %===Contains the segment indices on which the coupling conditions depend
Tj{k} = gseg{k}{2};
Deltaik = gseg{k}{3};
Deltai = [Deltai; Deltaik];

Gammaik = gseg{k}{4};
gammai = [gammai; Gammaik(:)];

if length(gseg{k})>=5 && Xj{k}(1)==0
    
    data.coupling{k}.func = gseg{k}{5};
    if is_empty_or_func(gseg{k}{6})
        data.coupling{k}.func_dt = gseg{k}{6};
        if is_empty_or_func(gseg{k}{7})
            data.coupling{k}.func_dp = gseg{k}{7};
            if is_empty_or_func(gseg{k}{8})
                data.coupling{k}.func_dtdt = gseg{k}{8};
                if is_empty_or_func(gseg{k}{9})
                    data.coupling{k}.func_dtdp = gseg{k}{9};
                    if is_empty_or_func(gseg{k}{10})
                        data.coupling{k}.func_dpdp = gseg{k}{10};
                    end
                end
            end
        end
    end
    
elseif length(gseg{k})>=5
    data.coupling{k}.func = gseg{k}{5};
    if is_empty_or_func(gseg{k}{6})
        data.coupling{k}.func_dp = gseg{k}{6};
        if is_empty_or_func(gseg{k}{7})
            data.coupling{k}.func_dpdp = gseg{k}{7};
        end
    end
    
    if length(gseg{k})>7
        data.coupling{k}.rows = gseg{k}{8};
        data.coupling{k}.cols = gseg{k}{9};
    else
        data.coupling{k}.rows = [];
        data.coupling{k}.cols = [];
    end
else
    data.coupling{k} = [];
end
end

data.Xj = Xj;
data.Tj = Tj;

data.pnames = {};
if iscellstr(str.peek('cell'))
  data.pnames = str.get('cell');
end
p0 = str.get;

fname = sprintf(oid,data.xi);
fbid = coco_get_id(fname, 'ddaecoll'); % Create toolbox instance identifier
fdata    = coco_get_func_data(prob, fbid, 'data');
data.ddaecoll_seg = fdata;

sol = [p0; Deltai; gammai];

data = coupling_init_data(data,p0);        % Build toolbox data

prob = coupling_construct_seg(prob, tbid, data,sol); % Append continuation problem

end

function flag = is_empty_or_func(x)
%IS_EMPTY_OR_FUNC   Check if input is empty or contains a function handle.
flag = isempty(x) || isa(x, 'function_handle');
end

