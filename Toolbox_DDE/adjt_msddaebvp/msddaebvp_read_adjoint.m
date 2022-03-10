function [sol, data] = msddaebvp_read_adjoint(oid, run, varargin)
%COLL_READ_ADJOINT   Read adjoint data from disk.
%
% [SOL DATA] = COLL_READ_ADJOINT(VARARGIN)
%
% VARARGIN = { [OID] RUN LAB [POINT] }
% Read adjoint data from solution data file of run RUN with solution label
% LAB.
%
% VARARGIN = { '' '' DATA }
% Construct initial adjoint solution and data structure.
%
% On input:
%
% OID : Optional object instance identifier (string, optional).
% RUN : Run identifier (string or cell-array of strings).
% LAB : Solution label (integer).
%
% DATA : Data structure.
%
% On output:
%
% SOL  : Adjoint solution structure.
% DATA : Adjoint data structure or unchanged copy of input argument DATA.
%
% In the first calling form, COLL_READ_ADJOINT reconstructs the adjoint
% solution and toolbox data structures from a saved solution and constructs
% restart information if solution is a bifurcation point. More
% specifically, denote with
%
%   'seg'     :  a branch of trajectory segments.
%
% and with 'BR(SP)' a special point detected along a branch of trajectory
% segments of type BR. The solution structure SOL will have the fields
%
%   SOL.L     :  Lagrange multipliers for trajectory segment.
%   SOL.TL    :  Differentials of Lagrange multipliers for trajectory
%                segment.
%   SOL.T0_L  :  Lagrange multiplier for initial time.
%   SOL.T0_TL :  Differential of Lagrange multiplier for initial time.
%
% and additional fields encoding an initial solution point as required by
% COLL_CONSTRUCT_ADJT. Depending on the types of the solution branch, the
% return value of SOL will have the following additional fields:
%
%   'seg(BP)' : For branch-points the fields SOL.TL0 and SOL.T0_TL0  will be
%               initialized to singular vectors transversal to SOL.TL and
%               sol.T0_TL.
%
% See also: COCO_READ_ADJOINT, ADJT_ISOL2COLL, ADJT_COLL2COLL,
% COLL_ADJT_INIT_DATA

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: ep_read_solution.m 2902 2015-10-09 18:06:32Z hdankowicz $

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




