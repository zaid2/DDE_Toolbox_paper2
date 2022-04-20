function [sol,data] = coupling_read_adjoint(oid, W, run, varargin)
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

if isempty(W)
    tbid = coco_get_id(oid, 'coupling'); % Create toolbox instance identifier
else
    fname = sprintf('coupling%d',W);
    tbid = coco_get_id(oid,fname);
end

[data, chart1, lidx1] = coco_read_adjoint(tbid, run, lab, 'data', ...
    'chart', 'lidx');

seg  = data.coupling_seg.ddaecoll_seg;
dim = seg.dim;

Ci = data.coupling_seg.Ci;
dim_xieks = data.coupling_seg.dim_xiek;
dim_xibks = data.coupling_seg.dim_xibk;

l0_coupling = chart1.x(1:seg.xbpdim,1);
tl_coupling = chart1.t(1:seg.xbpdim,1);

l0_coupling = reshape(l0_coupling,seg.xbp_shp);
tl_coupling = reshape(tl_coupling,seg.xbp_shp);

sol.l0_coupling = l0_coupling(:,seg.tbp_idx);
sol.tl_coupling = tl_coupling(:,seg.tbp_idx);

idx_bc = seg.xbpdim+1:seg.xbpdim+(Ci+1)+(Ci-1)*dim+dim_xieks+dim_xibks;

sol.l0_bc = chart1.x(idx_bc(:),1);
sol.tl_bc = chart1.t(idx_bc(:),1);

if ~isempty(data.coupling_seg.pnames)
    [chart2, lidx2] = coco_read_adjoint(coco_get_id(tbid, 'pars'), run, ...
        lab, 'chart', 'lidx');
    sol.pars_l0 = chart2.x;
    sol.pars_tl = chart2.t;
else
    lidx2 = [];
    sol.pars_l0 = [];
    sol.pars_tl = [];
end

sol.tbp = seg.tbp(seg.tbp_idx);


if ~isempty(pt) && strcmpi(pt,'BP')
    cdata   = coco_get_chart_data(chart1, 'lsol');
    tl0 = cdata.v(lidx1);
    
    tl0_coupling = tl0(1:seg.xbpdim,1);
    tl0_coupling = reshape(tl0_coupling,seg.xbp_shp);
    sol.tl0_coupling = tl0_coupling(:,seg.tbp_idx);
    sol.tl0_bc  = tl0(idx_bc(:),1);
    
    sol.pars_tl0 = cdata.v(lidx2);
    
else
    sol.tl0 = [];
    sol.tl0_coupling = [];
    sol.pars_tl0 = [];
end

end

function [sol, data] = read_sol_from(data)
% Construct solution structure from data.

% sol.l0 = [];
sol.l0 = [];
sol.l0_coupling = [];
sol.tl0 = [];
sol.pars_l0 = [];
sol.pars_tl0 = [];
end
