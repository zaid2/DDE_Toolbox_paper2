function [sol, data] = msddaebvp_read_adjoint(oid, run, varargin)
% Copyright (C) Zaid Ahsan

if isempty(oid) && isempty(run)
  [sol, data] = read_sol_from(varargin{:});
  return
end

% if nargin<3
%   [oid, run, lab] = coco_deal('', oid, run);
% else
%   lab = varargin{1};
% end

% info = coco_read_tb_info(oid, run, lab, 'ddaecoll');
% format      = info.format;
% branch_type = info.branch_type;

lab = varargin{1};
if numel(varargin)>=2
    pt = varargin{2};
else
    pt = [];
end


tbid = coco_get_id(oid, 'msbvp');


%=====boundary conditions
    [data, chart1, lidx1] = coco_read_adjoint(tbid, run, lab, 'data', ...
      'chart', 'lidx');
    sol.bc_l0  = chart1.x;
    sol.bc_tl0 = chart1.t;

%=====parameter gluing conditions    
    for i=2:data.nsegs
     gfid = coco_get_id(tbid, sprintf('shared%d', i-1));
    [chart2, lidx] = coco_read_adjoint(gfid, run, lab, ...
        'chart', 'lidx');
    lidx2{i-1} = lidx;
    sol.glue_l0{i-1}  = chart2.x;
    sol.glue_tl0{i-1} = chart2.t;
    end
    
%=====adjoint variables associated with the problem parameters     
    if ~isempty(data.pnames)
    pfid   = coco_get_id(tbid, 'pars');    
    [chart3, lidx3] = coco_read_adjoint(pfid, run, ...
        lab, 'chart', 'lidx');
    sol.pars_l0  = chart3.x;
    sol.pars_tl0 = chart3.t;
    else
        lidx3 = [];
        sol.pars_l0  = [];
        sol.pars_tl0 = [];
    end

    if ~isempty(pt) && strcmpi(pt,'BP')
        chart = coco_read_solution(run, lab, 'chart');
        cdata = coco_get_chart_data(chart, 'lsol');        
        sol.bc_tl0      = cdata.v(lidx1);
        for i=2:data.nsegs
           sol.glue_tl0{i-1} = cdata.v(lidx2{i-1}); 
        end
        sol.pars_tl0 = cdata.v(lidx3);
    end
    
%======Adjoints for the segments
sol.msbvp = cell(1,data.nsegs);
for i=1:data.nsegs
segoid = coco_get_id(tbid,sprintf('seg%d',i));
sol.msbvp{i} = ddaecoll_read_adjoint(segoid,run,lab);
end

end

function [sol, data] = read_sol_from(data)
% Construct solution structure from data.

sol.bc_l0       = [];
sol.bc_tl0      = [];
for i=2:data.nsegs
sol.glue_l0{i-1}  = [];
sol.glue_tl0{i-1} = [];
end
sol.pars_l0 = [];
sol.pars_tl0 = [];
sol.msbvp = [];
end




