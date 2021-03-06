function prob = msddaebvp_close_segs(prob, tbid, data)
% Copyright (C) Zaid Ahsan

if ~isempty(data.bc_update) % Optional inclusion of boundary conditions function data
  data.tbid = tbid;
  data = coco_func_data(data); % Convert to func_data class for shared access
  prob = coco_add_slot(prob, tbid, @msbvp_bc_update, data, 'update');
end
T0_idx = zeros(data.nsegs,1);
T_idx  = zeros(data.nsegs,1);
x0_idx = [];
x1_idx = [];
s_idx  = cell(1, data.nsegs);
for i=1:data.nsegs
  fid      = coco_get_id(tbid,sprintf('seg%d.ddaecoll', i)); % Create 'coll' toolbox instance identifier
  [fdata uidx] = coco_get_func_data(prob, fid, 'data', 'uidx'); % Extract 'coll' data structure and context-dependent index array
  T0_idx(i)= uidx(fdata.T0_idx);           % Subset to T0
  T_idx(i) = uidx(fdata.T_idx);            % Subset to T
  x0_idx   = [x0_idx; uidx(fdata.x0_idx)]; % Subset to x0
  x1_idx   = [x1_idx; uidx(fdata.x1_idx)]; % Subset to x1
  s_idx{i} = uidx(fdata.p_idx);            % Subset to p
end
uidx = [T0_idx; T_idx; x0_idx; x1_idx; s_idx{1}]; % Use only one copy of problem parameters
data.uidx = uidx;
if isempty(data.dfbcdxhan) % Optional inclusion of explicit Jacobian of boundary conditions
  prob = coco_add_func(prob, tbid, @msbvp_F, data, ...
    'zero', 'uidx', uidx);
else
  prob = coco_add_func(prob, tbid, @msbvp_F, @msbvp_DFDU, data, ...
    'zero', 'uidx', uidx);
end
for i=2:data.nsegs % Glue redundant copies of problem parameters
  fid  = coco_get_id(tbid, sprintf('shared%d', i-1));
  prob = coco_add_glue(prob, fid, s_idx{1}, s_idx{i});
end
if ~isempty(data.pnames) % Optional monitor functions
  fid  = coco_get_id(tbid, 'pars');
  prob = coco_add_pars(prob, fid, s_idx{1}, data.pnames);
end
prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');

end

function [data y] = msbvp_F(prob, data, u)
%MSBVP_F   COCO-compatible wrapper to boundary conditions.
%
% Expects encoding of boundary conditions as function of T, x0, x1, and p.

T0 = u(data.T0_idx); % Extract initial times
T  = u(data.T_idx);  % Extract interval lengths
x0 = u(data.x0_idx); % Extract trajectory end points at t=0
x1 = u(data.x1_idx); % Extract trajectory end points at t=1
p  = u(data.p_idx);  % Extract single copy of problem parameters

y  = data.fbchan(data.bc_data, T0, T, x0, x1, p);

end

function [data J] = msbvp_DFDU(prob, data, u)
%MSBVP_DFDU   COCO-compatible wrapper to linearization of boundary conditions.
%
% Expects encoding of boundary conditions as function of T, x0, x1, and p.

T0 = u(data.T0_idx); % Extract initial times
T  = u(data.T_idx);  % Extract interval lengths
x0 = u(data.x0_idx); % Extract trajectory end points at t=0
x1 = u(data.x1_idx); % Extract trajectory end points at t=1
p  = u(data.p_idx);  % Extract single copy of problem parameters

J  = data.dfbcdxhan(data.bc_data, T0, T, x0, x1, p);
% if p(3)>0.1
%     [data, J2] = coco_ezDFDX('f(o,d,x)', prob, data, @msbvp_F, u);
% end
% [data, dJ] = coco_ezDFDX('f(o,d,x)', prob, data, @msbvp_F, u);
end

function data = msbvp_bc_update(prob, data, cseg, varargin)
%MSBVP_BC_UPDATE   COCO-compatible wrapper to boundary condition function data update function.
%
% Use information about current solution to update function data
% parameterizing the execution of the boundary condition zero function.

uidx = coco_get_func_data(prob, data.tbid, 'uidx'); % Context-dependent index set
u    = cseg.src_chart.x(uidx); % Current chart
T0   = u(data.T0_idx); % Extract initial times
T    = u(data.T_idx);  % Extract interval lengths
x0   = u(data.x0_idx); % Extract trajectory end points at t=0
x1   = u(data.x1_idx); % Extract trajectory end points at t=1
p    = u(data.p_idx);  % Extract problem parameters
data.bc_data = data.bc_update(data.bc_data, T0, T, x0, x1, p);

end










