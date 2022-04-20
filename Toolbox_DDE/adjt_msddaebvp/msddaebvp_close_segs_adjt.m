function prob = msddaebvp_close_segs_adjt(prob, tbid, data, sol)
% Copyright (C) Zaid Ahsan

msbvp_data = data.msbvp_data;

T0_aidx = zeros(msbvp_data.nsegs,1);
T_aidx  = zeros(msbvp_data.nsegs,1);
x0_aidx = [];
x1_aidx = [];
% s_aidx  = cell(1, data.nsegs);
for i=1:msbvp_data.nsegs
  fid      = coco_get_id(tbid,sprintf('seg%d.ddaecoll', i)); % Create 'coll' toolbox instance identifier
  [adata, axidx] = coco_get_adjt_data(prob, fid, 'data', 'axidx'); % Extract 'coll' data structure and context-dependent index array
  opt = adata.ddaecoll_opt;
  T0_aidx(i)= axidx(opt.T0_idx);           % Subset to T0
  T_aidx(i) = axidx(opt.T_idx);            % Subset to T
  x0_aidx   = [x0_aidx; axidx(opt.x0_idx)]; % Subset to x0
  x1_aidx   = [x1_aidx; axidx(opt.x1_idx)]; % Subset to x1
  s_aidx{i} = axidx(opt.p_idx);            % Subset to p
end
aidx = [T0_aidx; T_aidx; x0_aidx; x1_aidx; s_aidx{1}]; % Use only one copy of problem parameters
data.aidx = aidx;

prob = coco_add_adjt(prob, tbid, 'aidx', aidx, 'l0',sol.bc_l0,'tl0',...
                                                        sol.bc_tl0); 
for i=2:msbvp_data.nsegs
   gfid = coco_get_id(tbid, sprintf('shared%d', i-1));
   prob = coco_add_adjt(prob, gfid, 'aidx', [s_aidx{1};s_aidx{i}], 'l0',...
                              sol.glue_l0{i-1},'tl0',sol.glue_tl0{i-1});
end

if ~isempty(msbvp_data.pnames) % Optional monitor functions
  pfid   = coco_get_id(tbid, 'pars');
  dnames = coco_get_id('d', msbvp_data.pnames);
  
%   prob = coco_add_pars(prob, fid, s_idx{1}, data.pnames);
  prob   = coco_add_adjt(prob, pfid, dnames, 'aidx', s_aidx{1}, ...
    'l0', sol.pars_l0, 'tl0', sol.pars_tl0);
end
% prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');

end

