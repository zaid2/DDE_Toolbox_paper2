% coco_use_recipes_toolbox ddecoll_vz1
                             
% Copyright (C) Zaid Ahsan                             
%% 2-dimensional continuation  
alpha = 1;
[t0,x0,y0] = init_guess_dde23(alpha);

prob = coco_prob();

prob = coco_set(prob, 'ddaecoll', 'NTST', 10);
prob = ddaecoll_isol2seg(prob,'seg1',@scalar_f,@scalar_dfdt,@scalar_dfdx,...
    @scalar_dfdy,@scalar_dfdp,@scalar_dfdtdt,@scalar_dfdtdx,...
    @scalar_dfdtdy,@scalar_dfdtdp,@scalar_dfdxdx,@scalar_dfdxdy,...
    @scalar_dfdxdp,@scalar_dfdydy,@scalar_dfdydp,@scalar_dfdpdp,...
    t0,x0,y0,{'p1' 'p2' 'p3' 'p4' 'p5' 'p6' 'p7' 'p8'},...
    [0;0;0;0;0;0;0;0]);
prob = adjt_ddaecoll_isol2seg(prob, 'seg1');

% Adding the boundary conditions
fbid = coco_get_id('seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata, uidx] = coco_get_func_data(prob, fbid, 'data', 'uidx');
prob = coco_add_func(prob,'bc',@bc,@bc_du, fdata, 'zero','uidx',...
    [uidx(fdata.T0_idx); uidx(fdata.T_idx)]);
[data, axidx] = coco_get_adjt_data(prob, 'seg1.ddaecoll', 'data', 'axidx');
opt  = data.ddaecoll_opt;
prob = coco_add_adjt(prob, 'bc', 'aidx', ...
   axidx([opt.T0_idx; opt.T_idx]));


% Coupling conditions

% 1. Coupling condition 1: Delay Term
g11 = {0,0,1,[0,1/2],@hist,@hist_dt,@hist_dp,@hist_dtdt,@hist_dtdp,@hist_dpdp};
g12 = {1,1,1/2,[1/2,1]};
xi = 1;
yiid = 1;
W = 1;
prob = coupling_isol2seg(prob,'seg1',W,xi,yiid,{g11,g12});  
prob = adjt_coupling_isol2seg(prob, 'seg1',W);

[data_g1, uidx_g1] = coco_get_func_data(prob, 'seg1.coupling1', 'data', 'uidx');

% Adding the corresp. gluing conditions
prob = coco_add_func(prob,'gluing1',@glue_coup1,@glue_coup1_du,[],'zero','uidx',...
            uidx_g1(data_g1.gamma_idx(2)));

[cdata1, cxidx1] = coco_get_adjt_data(prob, 'seg1.coupling1', 'data', 'axidx');
copt1 = cdata1.coupling_opt; 
axidx_glue = cxidx1(copt1.gamma_adjt_idx(2)); 
prob = coco_add_adjt(prob, 'gluing1', 'aidx', axidx_glue);         

% 2. Coupling condition 2: Controls
g2 = {0,0, 1,[0,1],@fcn_ctrl,@fcn_ctrl_dt,@fcn_ctrl_dp,@fcn_ctrl_dtdt,@fcn_ctrl_dtdp,@fcn_ctrl_dpdp};
xi = 1;
yiid = 2;
W = 2;
prob = coupling_isol2seg(prob,'seg1',W,xi,yiid,{g2});                                 
prob = adjt_coupling_isol2seg(prob, 'seg1',W);

[data_g2, uidx_g2] = coco_get_func_data(prob, 'seg1.coupling2', 'data', 'uidx');

% Adding the corresp. gluing conditions
prob = coco_add_func(prob,'gluing2',@glue_coup2,@glue_coup2_du,[],'zero','uidx',...
            uidx_g2(data_g2.Delta_idx(1)));
[cdata1, cxidx1] = coco_get_adjt_data(prob, 'seg1.coupling2', 'data', 'axidx');
copt1 = cdata1.coupling_opt; 
axidx_glue = cxidx1(copt1.Delta_adjt_idx); 
prob = coco_add_adjt(prob, 'gluing2', 'aidx', axidx_glue);         

% Objective function
fbid = coco_get_id('seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata, uidx] = coco_get_func_data(prob, fbid, 'data', 'uidx');
uxidx = [fdata.xbp_idx; fdata.ybp_idx(2:2:end)];
prob = coco_add_func(prob,'obj',@sys_obj,fdata,'inactive','obj','uidx', uxidx);

[data, axidx] = coco_get_adjt_data(prob, 'seg1.ddaecoll', 'data', 'axidx');
opt  = data.ddaecoll_opt;
aidx = [opt.xcn_idx; opt.ybp_idx(2:2:end)];
prob = coco_add_adjt(prob,'obj',@sys_obj_adj,fdata,'d.obj','aidx',aidx);

prob = coco_set(prob,'cont','atlas','kd','PtMX',40,...
    'NAdapt',0,'R',50,'R_max',100,'R_min',1e-4, 'MaxRes', 100,'al_max',15);

prob = coco_add_event(prob,'Event1','d.obj',1);

bd1 = coco(prob, 'run1', [], 2,{'obj' 'd.obj' 'p1' 'd.p1' 'd.p2' 'd.p3' 'd.p4' 'd.p5' 'd.p6' 'd.p7' 'd.p8'},...
                                 {[],[-1,1],[-2.5,0.5],[]});


%% Run to drive d.p1 to zero

run = 'run1';
labs = coco_bd_labs(bd1,'Event1');
lab   = labs(end);

% branch switch data
% chart = coco_read_solution(run, lab, 'chart');
% cdata = coco_get_chart_data(chart, 'lsol');

prob = coco_prob();
prob = ddaecoll_sol2seg(prob, 'seg1', run, lab); % Reconstruct 'coll' continuation problem
prob = adjt_ddaecoll_sol2seg(prob, 'seg1',run,lab);

% Adding the boundary conditions
fbid = coco_get_id('seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata, uidx] = coco_get_func_data(prob, fbid, 'data', 'uidx');
prob = coco_add_func(prob,'bc',@bc,@bc_du, fdata, 'zero','uidx',...
    [uidx(fdata.T0_idx); uidx(fdata.T_idx)]);
[data, axidx] = coco_get_adjt_data(prob, 'seg1.ddaecoll', 'data', 'axidx');
opt  = data.ddaecoll_opt;
[chart,aidx] = coco_read_adjoint('bc', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'bc', 'aidx', ...
   axidx([opt.T0_idx; opt.T_idx]),'l0',chart.x);

% Coupling condition 1:
W = 1;
prob = coupling_sol2seg(prob, 'seg1', W, run, lab);
prob = adjt_coupling_sol2seg(prob, 'seg1', W, run, lab);

[data_g1, uidx_g1] = coco_get_func_data(prob, 'seg1.coupling1', 'data', 'uidx');

% Adding the corresp. gluing conditions for coupling 1
prob = coco_add_func(prob,'gluing1',@glue_coup1,@glue_coup1_du,[],'zero','uidx',...
            uidx_g1(data_g1.gamma_idx(2)));

[cdata1, cxidx1] = coco_get_adjt_data(prob, 'seg1.coupling1', 'data', 'axidx');
copt1 = cdata1.coupling_opt; 
axidx_glue = cxidx1(copt1.gamma_adjt_idx(2)); 
[chart, aidx] = coco_read_adjoint('gluing1', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'gluing1', 'aidx', axidx_glue,'l0',chart.x); 

% Coupling condition 2:
W = 2;
prob = coupling_sol2seg(prob, 'seg1', W, run, lab);
prob = adjt_coupling_sol2seg(prob, 'seg1',W,run,lab);

% Adding the corresp. gluing conditions
[data_g2, uidx_g2] = coco_get_func_data(prob, 'seg1.coupling2', 'data', 'uidx');
prob = coco_add_func(prob,'gluing2',@glue_coup2,@glue_coup2_du,[],'zero','uidx',...
            uidx_g2(data_g2.Delta_idx(1)));
[cdata1, cxidx1] = coco_get_adjt_data(prob, 'seg1.coupling2', 'data', 'axidx');
copt1 = cdata1.coupling_opt; 
axidx_glue = cxidx1(copt1.Delta_adjt_idx); 
[chart, aidx] = coco_read_adjoint('gluing2', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'gluing2', 'aidx', axidx_glue,'l0',chart.x); 

% Objective function
fbid = coco_get_id('seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata, uidx] = coco_get_func_data(prob, fbid, 'data', 'uidx');
uxidx = [fdata.xbp_idx; fdata.ybp_idx(2:2:end)];
prob = coco_add_func(prob,'obj',@sys_obj,fdata,'inactive','obj','uidx', uxidx);

[data, axidx] = coco_get_adjt_data(prob, 'seg1.ddaecoll', 'data', 'axidx');
opt  = data.ddaecoll_opt;
aidx_obj = [opt.xcn_idx; opt.ybp_idx(2:2:end)];
[chart, aidx] = coco_read_adjoint('obj', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob,'obj',@sys_obj_adj,fdata,'d.obj','aidx',aidx_obj,'l0',chart.x);

axidx = coco_get_adjt_data( prob, 'obj', 'afidx');
fcn  = @(f) @(p,d,u,l,v) deal(d, f(u,l,v));
prob = coco_add_comp(prob, 'fun1', fcn(@(u,l,v) l-1), [], 'zero', ...
   'lidx', axidx);

axidx = coco_get_adjt_data( prob, 'seg1.ddaecoll.pars', 'afidx');
fcn  = @(f) @(p,d,u,l,v) deal(d, f(u,l,v));
prob = coco_add_comp(prob, 'fun2', fcn(@(u,l,v) l-0), [], 'zero', ...
   'lidx', axidx(1));


prob = coco_set(prob,'cont','ItMX',150);
prob = coco_set(prob,'corr','ItMX',30);

prob = coco_add_event(prob,'Event2','d.p2',1e-6);

bd2 = coco(prob, 'run2', [], 1,{'d.p2' 'p1' 'p2' 'obj' 'd.obj' 'd.p1' 'd.p3'...
                                       'd.p4' 'd.p5' 'd.p6' 'd.p7' 'd.p8'},{[0,5]});

%% run to drive d.p3 to zero

run = 'run2';
labs = coco_bd_labs(bd2,'Event2');
lab   = labs(end);


% branch switch data
% chart = coco_read_solution(run, lab, 'chart');
% cdata = coco_get_chart_data(chart, 'lsol');

prob = coco_prob();
prob = ddaecoll_sol2seg(prob, 'seg1', run, lab); % Reconstruct 'coll' continuation problem
prob = adjt_ddaecoll_sol2seg(prob, 'seg1',run,lab);

% Adding the boundary conditions
fbid = coco_get_id('seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata, uidx] = coco_get_func_data(prob, fbid, 'data', 'uidx');
prob = coco_add_func(prob,'bc',@bc,@bc_du, fdata, 'zero','uidx',...
    [uidx(fdata.T0_idx); uidx(fdata.T_idx)]);
[data, axidx] = coco_get_adjt_data(prob, 'seg1.ddaecoll', 'data', 'axidx');
opt  = data.ddaecoll_opt;
[chart,aidx] = coco_read_adjoint('bc', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'bc', 'aidx', ...
   axidx([opt.T0_idx; opt.T_idx]),'l0',chart.x);

% Coupling condition 1:
W = 1;
prob = coupling_sol2seg(prob, 'seg1', W, run, lab);
prob = adjt_coupling_sol2seg(prob, 'seg1', W, run, lab);

[data_g1, uidx_g1] = coco_get_func_data(prob, 'seg1.coupling1', 'data', 'uidx');

% Adding the corresp. gluing conditions for coupling 1
prob = coco_add_func(prob,'gluing1',@glue_coup1,@glue_coup1_du,[],'zero','uidx',...
            uidx_g1(data_g1.gamma_idx(2)));

[cdata1, cxidx1] = coco_get_adjt_data(prob, 'seg1.coupling1', 'data', 'axidx');
copt1 = cdata1.coupling_opt; 
axidx_glue = cxidx1(copt1.gamma_adjt_idx(2)); 
[chart, aidx] = coco_read_adjoint('gluing1', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'gluing1', 'aidx', axidx_glue,'l0',chart.x); 

% Coupling condition 2:
W = 2;
prob = coupling_sol2seg(prob, 'seg1', W, run, lab);
prob = adjt_coupling_sol2seg(prob, 'seg1',W,run,lab);

% Adding the corresp. gluing conditions
[data_g2, uidx_g2] = coco_get_func_data(prob, 'seg1.coupling2', 'data', 'uidx');
prob = coco_add_func(prob,'gluing2',@glue_coup2,@glue_coup2_du,[],'zero','uidx',...
            uidx_g2(data_g2.Delta_idx(1)));
[cdata1, cxidx1] = coco_get_adjt_data(prob, 'seg1.coupling2', 'data', 'axidx');
copt1 = cdata1.coupling_opt; 
axidx_glue = cxidx1(copt1.Delta_adjt_idx); 
[chart, aidx] = coco_read_adjoint('gluing2', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'gluing2', 'aidx', axidx_glue,'l0',chart.x); 

% Objective function
fbid = coco_get_id('seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata, uidx] = coco_get_func_data(prob, fbid, 'data', 'uidx');
uxidx = [fdata.xbp_idx; fdata.ybp_idx(2:2:end)];
prob = coco_add_func(prob,'obj',@sys_obj,fdata,'inactive','obj','uidx', uxidx);

[data, axidx] = coco_get_adjt_data(prob, 'seg1.ddaecoll', 'data', 'axidx');
opt  = data.ddaecoll_opt;
aidx_obj = [opt.xcn_idx; opt.ybp_idx(2:2:end)];
[chart, aidx] = coco_read_adjoint('obj', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob,'obj',@sys_obj_adj,fdata,'d.obj','aidx',aidx_obj,'l0',chart.x);

prob = coco_set(prob,'cont','ItMX',150);
bd3 = coco(prob, 'run3', [], 1,{'d.p3' 'obj' 'p1' 'p2' 'p3'...
                                       'd.p4' 'd.p5' 'd.p6' 'd.p7' 'd.p8'},{[-0.4,0]});

%% run to drive d.p4 to zero

run = 'run3';
lab = 4;

prob = coco_prob();
prob = ddaecoll_sol2seg(prob, 'seg1', run, lab); % Reconstruct 'coll' continuation problem
prob = adjt_ddaecoll_sol2seg(prob, 'seg1',run,lab);

% Adding the boundary conditions
fbid = coco_get_id('seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata, uidx] = coco_get_func_data(prob, fbid, 'data', 'uidx');
prob = coco_add_func(prob,'bc',@bc,@bc_du, fdata, 'zero','uidx',...
    [uidx(fdata.T0_idx); uidx(fdata.T_idx)]);
[data, axidx] = coco_get_adjt_data(prob, 'seg1.ddaecoll', 'data', 'axidx');
opt  = data.ddaecoll_opt;
[chart,aidx] = coco_read_adjoint('bc', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'bc', 'aidx', ...
   axidx([opt.T0_idx; opt.T_idx]),'l0',chart.x);

% Coupling condition 1:
W = 1;
prob = coupling_sol2seg(prob, 'seg1', W, run, lab);
prob = adjt_coupling_sol2seg(prob, 'seg1', W, run, lab);

[data_g1, uidx_g1] = coco_get_func_data(prob, 'seg1.coupling1', 'data', 'uidx');

% Adding the corresp. gluing conditions for coupling 1
prob = coco_add_func(prob,'gluing1',@glue_coup1,@glue_coup1_du,[],'zero','uidx',...
            uidx_g1(data_g1.gamma_idx(2)));

[cdata1, cxidx1] = coco_get_adjt_data(prob, 'seg1.coupling1', 'data', 'axidx');
copt1 = cdata1.coupling_opt; 
axidx_glue = cxidx1(copt1.gamma_adjt_idx(2)); 
[chart, aidx] = coco_read_adjoint('gluing1', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'gluing1', 'aidx', axidx_glue,'l0',chart.x); 

% Coupling condition 2:
W = 2;
prob = coupling_sol2seg(prob, 'seg1', W, run, lab);
prob = adjt_coupling_sol2seg(prob, 'seg1',W,run,lab);

% Adding the corresp. gluing conditions
[data_g2, uidx_g2] = coco_get_func_data(prob, 'seg1.coupling2', 'data', 'uidx');
prob = coco_add_func(prob,'gluing2',@glue_coup2,@glue_coup2_du,[],'zero','uidx',...
            uidx_g2(data_g2.Delta_idx(1)));
[cdata1, cxidx1] = coco_get_adjt_data(prob, 'seg1.coupling2', 'data', 'axidx');
copt1 = cdata1.coupling_opt; 
axidx_glue = cxidx1(copt1.Delta_adjt_idx); 
[chart, aidx] = coco_read_adjoint('gluing2', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'gluing2', 'aidx', axidx_glue,'l0',chart.x); 

% Objective function
fbid = coco_get_id('seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata, uidx] = coco_get_func_data(prob, fbid, 'data', 'uidx');
uxidx = [fdata.xbp_idx; fdata.ybp_idx(2:2:end)];
prob = coco_add_func(prob,'obj',@sys_obj,fdata,'inactive','obj','uidx', uxidx);

[data, axidx] = coco_get_adjt_data(prob, 'seg1.ddaecoll', 'data', 'axidx');
opt  = data.ddaecoll_opt;
aidx_obj = [opt.xcn_idx; opt.ybp_idx(2:2:end)];
[chart, aidx] = coco_read_adjoint('obj', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob,'obj',@sys_obj_adj,fdata,'d.obj','aidx',aidx_obj,'l0',chart.x);

prob = coco_set(prob,'cont','ItMX',150);
bd4 = coco(prob, 'run4', [], 1,{'d.p4' 'obj' 'p1' 'p2' 'p3'...
                                       'p4' 'd.p5' 'd.p6' 'd.p7' 'd.p8'},{[0,8e-3]});

%% run to drive d.p5 to zero

run = 'run4';
lab   = 2;

% branch switch data
chart = coco_read_solution(run, lab, 'chart');
% cdata = coco_get_chart_data(chart, 'lsol');

prob = coco_prob();
prob = ddaecoll_sol2seg(prob, 'seg1', run, lab); % Reconstruct 'coll' continuation problem
prob = adjt_ddaecoll_sol2seg(prob, 'seg1',run,lab);

% Adding the boundary conditions
fbid = coco_get_id('seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata, uidx] = coco_get_func_data(prob, fbid, 'data', 'uidx');
prob = coco_add_func(prob,'bc',@bc,@bc_du, fdata, 'zero','uidx',...
    [uidx(fdata.T0_idx); uidx(fdata.T_idx)]);
[data, axidx] = coco_get_adjt_data(prob, 'seg1.ddaecoll', 'data', 'axidx');
opt  = data.ddaecoll_opt;
[chart,aidx] = coco_read_adjoint('bc', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'bc', 'aidx', ...
   axidx([opt.T0_idx; opt.T_idx]),'l0',chart.x);

% Coupling condition 1:
W = 1;
prob = coupling_sol2seg(prob, 'seg1', W, run, lab);
prob = adjt_coupling_sol2seg(prob, 'seg1', W, run, lab);

[data_g1, uidx_g1] = coco_get_func_data(prob, 'seg1.coupling1', 'data', 'uidx');

% Adding the corresp. gluing conditions for coupling 1
prob = coco_add_func(prob,'gluing1',@glue_coup1,@glue_coup1_du,[],'zero','uidx',...
            uidx_g1(data_g1.gamma_idx(2)));

[cdata1, cxidx1] = coco_get_adjt_data(prob, 'seg1.coupling1', 'data', 'axidx');
copt1 = cdata1.coupling_opt; 
axidx_glue = cxidx1(copt1.gamma_adjt_idx(2)); 
[chart, aidx] = coco_read_adjoint('gluing1', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'gluing1', 'aidx', axidx_glue,'l0',chart.x); 

% Coupling condition 2:
W = 2;
prob = coupling_sol2seg(prob, 'seg1', W, run, lab);
prob = adjt_coupling_sol2seg(prob, 'seg1',W,run,lab);

% Adding the corresp. gluing conditions
[data_g2, uidx_g2] = coco_get_func_data(prob, 'seg1.coupling2', 'data', 'uidx');
prob = coco_add_func(prob,'gluing2',@glue_coup2,@glue_coup2_du,[],'zero','uidx',...
            uidx_g2(data_g2.Delta_idx(1)));
[cdata1, cxidx1] = coco_get_adjt_data(prob, 'seg1.coupling2', 'data', 'axidx');
copt1 = cdata1.coupling_opt; 
axidx_glue = cxidx1(copt1.Delta_adjt_idx); 
[chart, aidx] = coco_read_adjoint('gluing2', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'gluing2', 'aidx', axidx_glue,'l0',chart.x); 

% Objective function
fbid = coco_get_id('seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata, uidx] = coco_get_func_data(prob, fbid, 'data', 'uidx');
uxidx = [fdata.xbp_idx; fdata.ybp_idx(2:2:end)];
prob = coco_add_func(prob,'obj',@sys_obj,fdata,'inactive','obj','uidx', uxidx);

[data, axidx] = coco_get_adjt_data(prob, 'seg1.ddaecoll', 'data', 'axidx');
opt  = data.ddaecoll_opt;
aidx_obj = [opt.xcn_idx; opt.ybp_idx(2:2:end)];
[chart, aidx] = coco_read_adjoint('obj', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob,'obj',@sys_obj_adj,fdata,'d.obj','aidx',aidx_obj,'l0',chart.x);



prob = coco_set(prob,'cont','ItMX',150);
bd5 = coco(prob, 'run5', [], 1,{'d.p5' 'obj' 'p1' 'p2' 'p3' 'p4'...
                                       'p5' 'd.p6' 'd.p7' 'd.p8'},{[0,2e-3]});


%% run to drive d.p6 to zero

run = 'run5';
lab   = 2;

prob = coco_prob();
prob = ddaecoll_sol2seg(prob, 'seg1', run, lab); % Reconstruct 'coll' continuation problem
prob = adjt_ddaecoll_sol2seg(prob, 'seg1',run,lab);

% Adding the boundary conditions
fbid = coco_get_id('seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata, uidx] = coco_get_func_data(prob, fbid, 'data', 'uidx');
prob = coco_add_func(prob,'bc',@bc,@bc_du, fdata, 'zero','uidx',...
    [uidx(fdata.T0_idx); uidx(fdata.T_idx)]);
[data, axidx] = coco_get_adjt_data(prob, 'seg1.ddaecoll', 'data', 'axidx');
opt  = data.ddaecoll_opt;
[chart,aidx] = coco_read_adjoint('bc', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'bc', 'aidx', ...
   axidx([opt.T0_idx; opt.T_idx]),'l0',chart.x);

% Coupling condition 1:
W = 1;
prob = coupling_sol2seg(prob, 'seg1', W, run, lab);
prob = adjt_coupling_sol2seg(prob, 'seg1', W, run, lab);

[data_g1, uidx_g1] = coco_get_func_data(prob, 'seg1.coupling1', 'data', 'uidx');

% Adding the corresp. gluing conditions for coupling 1
prob = coco_add_func(prob,'gluing1',@glue_coup1,@glue_coup1_du,[],'zero','uidx',...
            uidx_g1(data_g1.gamma_idx(2)));

[cdata1, cxidx1] = coco_get_adjt_data(prob, 'seg1.coupling1', 'data', 'axidx');
copt1 = cdata1.coupling_opt; 
axidx_glue = cxidx1(copt1.gamma_adjt_idx(2)); 
[chart, aidx] = coco_read_adjoint('gluing1', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'gluing1', 'aidx', axidx_glue,'l0',chart.x); 

% Coupling condition 2:
W = 2;
prob = coupling_sol2seg(prob, 'seg1', W, run, lab);
prob = adjt_coupling_sol2seg(prob, 'seg1',W,run,lab);

% Adding the corresp. gluing conditions
[data_g2, uidx_g2] = coco_get_func_data(prob, 'seg1.coupling2', 'data', 'uidx');
prob = coco_add_func(prob,'gluing2',@glue_coup2,@glue_coup2_du,[],'zero','uidx',...
            uidx_g2(data_g2.Delta_idx(1)));
[cdata1, cxidx1] = coco_get_adjt_data(prob, 'seg1.coupling2', 'data', 'axidx');
copt1 = cdata1.coupling_opt; 
axidx_glue = cxidx1(copt1.Delta_adjt_idx); 
[chart, aidx] = coco_read_adjoint('gluing2', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'gluing2', 'aidx', axidx_glue,'l0',chart.x); 

% Objective function
fbid = coco_get_id('seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata, uidx] = coco_get_func_data(prob, fbid, 'data', 'uidx');
uxidx = [fdata.xbp_idx; fdata.ybp_idx(2:2:end)];
prob = coco_add_func(prob,'obj',@sys_obj,fdata,'inactive','obj','uidx', uxidx);

[data, axidx] = coco_get_adjt_data(prob, 'seg1.ddaecoll', 'data', 'axidx');
opt  = data.ddaecoll_opt;
aidx_obj = [opt.xcn_idx; opt.ybp_idx(2:2:end)];
[chart, aidx] = coco_read_adjoint('obj', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob,'obj',@sys_obj_adj,fdata,'d.obj','aidx',aidx_obj,'l0',chart.x);



prob = coco_set(prob,'cont','ItMX',150);
bd6 = coco(prob, 'run6', [], 1,{'d.p6' 'obj' 'p1' 'p2' 'p3'...
                                       'p4' 'p5' 'p6' 'd.p7' 'd.p8'},{[-2.2e-3,0]});

%% run to drive d.p7 to zero

run = 'run6';
lab   = 4;

% branch switch data
chart = coco_read_solution(run, lab, 'chart');
% cdata = coco_get_chart_data(chart, 'lsol');

prob = coco_prob();
prob = ddaecoll_sol2seg(prob, 'seg1', run, lab); % Reconstruct 'coll' continuation problem
prob = adjt_ddaecoll_sol2seg(prob, 'seg1',run,lab);

% Adding the boundary conditions
fbid = coco_get_id('seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata, uidx] = coco_get_func_data(prob, fbid, 'data', 'uidx');
prob = coco_add_func(prob,'bc',@bc,@bc_du, fdata, 'zero','uidx',...
    [uidx(fdata.T0_idx); uidx(fdata.T_idx)]);
[data, axidx] = coco_get_adjt_data(prob, 'seg1.ddaecoll', 'data', 'axidx');
opt  = data.ddaecoll_opt;
[chart,aidx] = coco_read_adjoint('bc', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'bc', 'aidx', ...
   axidx([opt.T0_idx; opt.T_idx]),'l0',chart.x);

% Coupling condition 1:
W = 1;
prob = coupling_sol2seg(prob, 'seg1', W, run, lab);
prob = adjt_coupling_sol2seg(prob, 'seg1', W, run, lab);

[data_g1, uidx_g1] = coco_get_func_data(prob, 'seg1.coupling1', 'data', 'uidx');

% Adding the corresp. gluing conditions for coupling 1
prob = coco_add_func(prob,'gluing1',@glue_coup1,@glue_coup1_du,[],'zero','uidx',...
            uidx_g1(data_g1.gamma_idx(2)));

[cdata1, cxidx1] = coco_get_adjt_data(prob, 'seg1.coupling1', 'data', 'axidx');
copt1 = cdata1.coupling_opt; 
axidx_glue = cxidx1(copt1.gamma_adjt_idx(2)); 
[chart, aidx] = coco_read_adjoint('gluing1', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'gluing1', 'aidx', axidx_glue,'l0',chart.x); 

% Coupling condition 2:
W = 2;
prob = coupling_sol2seg(prob, 'seg1', W, run, lab);
prob = adjt_coupling_sol2seg(prob, 'seg1',W,run,lab);

% Adding the corresp. gluing conditions
[data_g2, uidx_g2] = coco_get_func_data(prob, 'seg1.coupling2', 'data', 'uidx');
prob = coco_add_func(prob,'gluing2',@glue_coup2,@glue_coup2_du,[],'zero','uidx',...
            uidx_g2(data_g2.Delta_idx(1)));
[cdata1, cxidx1] = coco_get_adjt_data(prob, 'seg1.coupling2', 'data', 'axidx');
copt1 = cdata1.coupling_opt; 
axidx_glue = cxidx1(copt1.Delta_adjt_idx); 
[chart, aidx] = coco_read_adjoint('gluing2', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'gluing2', 'aidx', axidx_glue,'l0',chart.x); 

% Objective function
fbid = coco_get_id('seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata, uidx] = coco_get_func_data(prob, fbid, 'data', 'uidx');
uxidx = [fdata.xbp_idx; fdata.ybp_idx(2:2:end)];
prob = coco_add_func(prob,'obj',@sys_obj,fdata,'inactive','obj','uidx', uxidx);

[data, axidx] = coco_get_adjt_data(prob, 'seg1.ddaecoll', 'data', 'axidx');
opt  = data.ddaecoll_opt;
aidx_obj = [opt.xcn_idx; opt.ybp_idx(2:2:end)];
[chart, aidx] = coco_read_adjoint('obj', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob,'obj',@sys_obj_adj,fdata,'d.obj','aidx',aidx_obj,'l0',chart.x);



prob = coco_set(prob,'cont','ItMX',150);
bd7 = coco(prob, 'run7', [], 1,{'d.p7' 'obj' 'p1' 'p2' 'p3'...
                                       'p4' 'p5' 'p6' 'p7' 'd.p8'},{[0,2e-3]});


%% run to drive d.p8 to zero

run = 'run7';
lab   = 2;

% branch switch data
chart = coco_read_solution(run, lab, 'chart');
% cdata = coco_get_chart_data(chart, 'lsol');

prob = coco_prob();
prob = ddaecoll_sol2seg(prob, 'seg1', run, lab); % Reconstruct 'coll' continuation problem
prob = adjt_ddaecoll_sol2seg(prob, 'seg1',run,lab);

% Adding the boundary conditions
fbid = coco_get_id('seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata, uidx] = coco_get_func_data(prob, fbid, 'data', 'uidx');
prob = coco_add_func(prob,'bc',@bc,@bc_du, fdata, 'zero','uidx',...
    [uidx(fdata.T0_idx); uidx(fdata.T_idx)]);
[data, axidx] = coco_get_adjt_data(prob, 'seg1.ddaecoll', 'data', 'axidx');
opt  = data.ddaecoll_opt;
[chart,aidx] = coco_read_adjoint('bc', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'bc', 'aidx', ...
   axidx([opt.T0_idx; opt.T_idx]),'l0',chart.x);

% Coupling condition 1:
W = 1;
prob = coupling_sol2seg(prob, 'seg1', W, run, lab);
prob = adjt_coupling_sol2seg(prob, 'seg1', W, run, lab);

[data_g1, uidx_g1] = coco_get_func_data(prob, 'seg1.coupling1', 'data', 'uidx');

% Adding the corresp. gluing conditions for coupling 1
prob = coco_add_func(prob,'gluing1',@glue_coup1,@glue_coup1_du,[],'zero','uidx',...
            uidx_g1(data_g1.gamma_idx(2)));

[cdata1, cxidx1] = coco_get_adjt_data(prob, 'seg1.coupling1', 'data', 'axidx');
copt1 = cdata1.coupling_opt; 
axidx_glue = cxidx1(copt1.gamma_adjt_idx(2)); 
[chart, aidx] = coco_read_adjoint('gluing1', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'gluing1', 'aidx', axidx_glue,'l0',chart.x); 

% Coupling condition 2:
W = 2;
prob = coupling_sol2seg(prob, 'seg1', W, run, lab);
prob = adjt_coupling_sol2seg(prob, 'seg1',W,run,lab);

% Adding the corresp. gluing conditions
[data_g2, uidx_g2] = coco_get_func_data(prob, 'seg1.coupling2', 'data', 'uidx');
prob = coco_add_func(prob,'gluing2',@glue_coup2,@glue_coup2_du,[],'zero','uidx',...
            uidx_g2(data_g2.Delta_idx(1)));
[cdata1, cxidx1] = coco_get_adjt_data(prob, 'seg1.coupling2', 'data', 'axidx');
copt1 = cdata1.coupling_opt; 
axidx_glue = cxidx1(copt1.Delta_adjt_idx); 
[chart, aidx] = coco_read_adjoint('gluing2', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'gluing2', 'aidx', axidx_glue,'l0',chart.x); 

% Objective function
fbid = coco_get_id('seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata, uidx] = coco_get_func_data(prob, fbid, 'data', 'uidx');
uxidx = [fdata.xbp_idx; fdata.ybp_idx(2:2:end)];
prob = coco_add_func(prob,'obj',@sys_obj,fdata,'inactive','obj','uidx', uxidx);

[data, axidx] = coco_get_adjt_data(prob, 'seg1.ddaecoll', 'data', 'axidx');
opt  = data.ddaecoll_opt;
aidx_obj = [opt.xcn_idx; opt.ybp_idx(2:2:end)];
[chart, aidx] = coco_read_adjoint('obj', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob,'obj',@sys_obj_adj,fdata,'d.obj','aidx',aidx_obj,'l0',chart.x);



prob = coco_set(prob,'cont','ItMX',150);
bd = coco(prob, 'run8', [], 1,{'d.p8' 'obj' 'p1' 'p2' 'p3'...
                                       'p4' 'p5' 'p6' 'p7' 'p8'},{[0,2e-3]});
                                   

%% Plotting the solution: 1

[sol,data] = ddaecoll_read_solution('seg1','run8',2);

%==plotting the state variables
figure(1)
plot(sol.t/2,sol.x,'Linewidth',1.5)

set(gca,'FontSize',12)
xlabel('$\tau$', 'Interpreter', 'Latex', 'Fontsize',14)
ylabel('$x$', 'Interpreter', 'Latex', 'Fontsize',14)
grid('on');
%==plotting the y variables
figure(2)
plot(sol.t/2,sol.y(:,1),sol.t/2,sol.y(:,2),'r--','Linewidth',1.5)
grid('on');
set(gca,'FontSize',12)
xlabel('$\tau$', 'Interpreter', 'Latex', 'Fontsize',14)
ylabel('$y(\tau)$', 'Interpreter', 'Latex', 'Fontsize',14)
legend('$y^{(1)}(\tau)$', '$y^{(2)}(\tau)$','Interpreter','Latex','Fontsize',14,'Edgecolor','None')

%% Plotting the solution: run 8

%===Lagrange multipliers and optimum solution
run = 'run8'; lab = 2;
[sol,data] = ddaecoll_read_solution('seg1',run,lab);

%==plotting the state variables
figure(1)
plot(sol.t/2,sol.x,'Linewidth',1.5)
grid('on')
set(gca,'FontSize',12)
xlabel('$\tau$', 'Interpreter', 'Latex', 'Fontsize',18)
ylabel('$x$', 'Interpreter', 'Latex', 'Fontsize',18)
title('State Variable')

%==plotting the y variables
figure(2)
plot(sol.t/2,sol.y(:,1),sol.t/2,sol.y(:,2),'r--','Linewidth',1.5)
grid('on')
set(gca,'FontSize',12)
xlabel('$\tau$', 'Interpreter', 'Latex', 'Fontsize',18)
ylabel('$y(\tau)$', 'Interpreter', 'Latex', 'Fontsize',18)
legend('$y^{(1)}(\tau)$', '$y^{(2)}(\tau)$','Interpreter','Latex','Fontsize',14,'Edgecolor','None')
title('Coupling Variables y')

%==plotting the Lagrange multipliers
figure(3)
sol = ddaecoll_read_adjoint('seg1',run,lab);
plot(sol.tbp,sol.l0,'Linewidth',1.5)
grid('on')
set(gca,'FontSize',12)
xlabel('$\tau$', 'Interpreter', 'Latex', 'Fontsize',18)
ylabel('$\lambda(\tau)$', 'Interpreter', 'Latex', 'Fontsize',18)
title('Lagrange Multiplier \lambda')

figure(4)
sol = coupling_read_adjoint('seg1',1,run,lab);
plot(sol.tbp,sol.l0_coupling,'Linewidth',1.5)
hold on
sol = coupling_read_adjoint('seg1',2,run,lab);
plot(sol.tbp,sol.l0_coupling,'r--','Linewidth',1.5)
grid('on');
set(gca,'FontSize',12)
xlabel('$\tau$', 'Interpreter', 'Latex', 'Fontsize',18)
ylabel('$\mu$', 'Interpreter', 'Latex', 'Fontsize',18)
legend('$\mu^{(1)}(\tau)$','$\mu^{(2)}(\tau)$','Interpreter','Latex','Fontsize',14,'Edgecolor','None')
title('Lagrange Multiplier \mu')

H=gca;
H.LineWidth=1.5;

% %====
% atlas = coco_bd_read('run_surface', 'atlas');
% hold on
% plot_atlas_kd(atlas.charts,3,2,4)
% grid('on')
% box('on')
% % ylim([0,1])
% 
% bd = coco_bd_read('run_trivial_etap1');
% col = coco_bd_col(bd,{'p1','d.obj','obj'});
% hold on
% plot3(col(1,:)',col(2,:)',col(3,:)','b--','LineWidth',2)

% bd = coco_bd_read('run_nontrivial_etap1');
% col = coco_bd_col(bd,{'p1','d.obj','obj'});
% hold on
% plot3(col(1,:)',col(2,:)',col(3,:)','r-','LineWidth',3)

