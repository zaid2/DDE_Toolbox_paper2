function sol = ddaecoll_read_adjoint(oid, run, varargin)
% Copyright (C) Zaid Ahsan, Mingwu Li

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


tbid = coco_get_id(oid, 'ddaecoll');

% switch format
%   
%   case 'ddaecoll.v1'
    
    [data, chart1, lidx1] = coco_read_adjoint(tbid, run, lab, 'data', ...
      'chart', 'lidx');
    
    seg  = data.ddaecoll_seg;
    
    sol.tbp = seg.tbp(seg.tbp_idx);
    l0      = reshape(chart1.x, seg.xbp_shp);
    sol.l0  = l0(:, seg.tbp_idx);
    tl      = reshape(chart1.t, seg.xbp_shp);
    sol.tl  = tl(:, seg.tbp_idx);
    sol.tl0 = [];
    
    if ~isempty(seg.pnames)
      [chart2, lidx2] = coco_read_adjoint(coco_get_id(tbid, 'pars'), run, ...
        lab, 'chart', 'lidx');
      sol.pars_l0 = chart2.x;
      sol.pars_tl = chart2.t;
    else
      lidx2 = [];
      sol.pars_l0 = [];
      sol.pars_tl = [];
    end
    sol.pars_tl0 = [];
    
    
    if ~isempty(pt) && strcmpi(pt,'BP')
        cdata = coco_get_chart_data(chart1, 'lsol');
        tl0      = reshape(cdata.v(lidx1), seg.xbp_shp);
        sol.tl0  = tl0(:, seg.tbp_idx);
        sol.pars_tl0 = cdata.v(lidx2);            
    end
    

end

function [sol, data] = read_sol_from(data)
% Construct solution structure from data.

sol.l0       = [];
sol.tl0      = [];
sol.pars_l0  = [];
sol.pars_tl0 = [];
% sol.T0_l0    = [];
% sol.T0_tl0   = [];

end
