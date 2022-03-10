% coco_use_recipes_toolbox ddecoll_vz1

%% Adding the vector field

 
% Construct 'ddecoll' arguments
alpha = 0; 
om2 = 1.5;  %===om
om1 = 1;    %===Om
Texc = 2*pi/om2;
varrho = 1/1.51111;
N = 5;
vphi = 2*pi*linspace(0,1,2*N+2);  %===== discretization in the phase

tau  = 2*pi/om2*linspace(0,1,10*(2*N+1))';   %===preparing the initial guess
rho1  = (1+om2^2)./(1+om2^2-cos(om2*tau)-om2*sin(om2*tau));
rho2  = (1+om2^2)./(1+om2^2-cos(om2*(tau))-om2*sin(om2*(tau)));

prob = coco_prob();
prob = coco_set(prob,'ddaecoll','NTST',10);
%=====Constructing the 2N+1 segments
ddaecoll_args = {};
x0 = [];
for i=1:2*N+1
  upx = repmat(rho1, [1 2]).*[cos(om1*tau+vphi(i)) sin(om1*tau+vphi(i))];
  upy = repmat(rho2, [1 2]).*[cos(om1*(tau)+vphi(i)) sin(om1*(tau)+vphi(i))];
  ddaecoll_args_temp = {@ddetorus @ddetorus_dfdt @ddetorus_dfdx @ddetorus_dfdy...
                               @ddetorus_dfdp @ddetorus_dfdtdt @ddetorus_dfdtdx...
                               @ddetorus_dfdtdy @ddetorus_dfdtdp @ddetorus_dfdxdx @ddetorus_dfdxdy...
                               @ddetorus_dfdxdp @ddetorus_dfdydy @ddetorus_dfdydp @ddetorus_dfdpdp...
                               tau upx upy [om1;Texc;alpha]};
  x0 = [x0; upx(1,:)'];                         
  ddaecoll_args = [ddaecoll_args ddaecoll_args_temp];
end

% Construct boundary conditions data, Fourier transform and rotation matrix
Th = 2*pi*(0:2*N)/(2*N+1);
Th = kron(1:N, Th');
F  = [ones(2*N+1,1) 2*reshape([cos(Th);sin(Th)], [2*N+1 2*N])]'/(2*N+1);

varrho = 1/1.51111;
Th  = (1:N)*2*pi*varrho;
SIN = [ zeros(size(Th)) ; sin(Th) ];
R   = diag([1 kron(cos(Th), [1 1])]);
R   = R  + diag(SIN(:), +1)- diag(SIN(:), -1);

Th  = -(1:N)*2 *pi*varrho;     %===difference between R and R1 is R1 has -2*pi*varrho
SIN = [ zeros(size(Th)) ; sin(Th) ];
R1   = diag([1 kron(cos(Th), [1 1])]);
R1   = R1  + diag(SIN(:), +1)- diag(SIN(:), -1);

data    = struct();
data.F  = kron(F, eye(2));         %====storing the Fourier discretization matrix
data.RF = kron(R*F, eye(2));       %====storing the Rotation matrix
data.R1F = kron(R1*F, eye(2));       %====storing the Rotation matrix
data.nsegs = 2*N+1;
data.dim = 2;

data = Update_x0(data, [], [], x0, [], []);  %===for phase condition update
prob = msddaebvp_isol2segs(prob, '', ddaecoll_args{:}, {'om1' 'Texc' 'alpha'},...
              @bc_ddetorus, @bc_ddetorus_du, data, @Update_x0); % Build 'coll' continuation problem

T = 2*pi/om2;
dim = 2;
uidx_glue = [];
%===Coupling Conditions
for i=1:2*N+1
xi = i;
g11 = {(1:2*N+1)',i,-(1-alpha/T), [0,alpha/T], @Mat_coupling,@Mat_coupling_dp,@Mat_coupling_dpdp,...
    1+(i-1)*dim:2*i,1:dim*(2*N+1)}; 
% g11 = [j_indices, Delta, [gamma_ikb,gamma_ike], Matrix, [],[],rows,cols]. Empty brackets are for derivatives of Matrix w.r.t p
g12 = {i, i, (alpha/T), [alpha/T,1]};

fid = sprintf('msbvp.seg%d', i);
prob = coupling_isol2seg(prob, fid, [], xi, 1:2, {g11,g12});

fbid = sprintf('msbvp.seg%d.coupling', i);
[data_gi, uidx_gi] = coco_get_func_data(prob, fbid, 'data', 'uidx');

uidx_glue = [uidx_glue; uidx_gi(data_gi.gamma_idx(2))];

end

fbid = coco_get_id('msbvp.seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata1, uidx1] = coco_get_func_data(prob, fbid, 'data', 'uidx');

uidx_glue = [uidx_glue; uidx1(fdata1.p_idx(3)); uidx1(fdata1.T_idx)];

prob = coco_add_func(prob,'gluing',@glue_coup,@glue_coup_du,data,'zero',...
                   'uidx', uidx_glue);

prob = coco_set(prob,'cont','atlas','kd','PtMX',40000,...
    'NAdapt',0,'R',.01,'R_max',250,'R_min',1e-4, 'MaxRes', 10,'al_max',10);

prob = coco_set(prob,'corr','ItMX',100);
% coco(prob,'run_test1',[],1,{'om1','om2','alpha'},{[0.8,6.5]})
coco(prob,'run_test0',[],1,{'alpha','om1','Texc'},{[0,0.1],[]})


%% ======= Run to drive the delay to unity
run = 'run_test0'; lab = 2;

prob = coco_prob();

prob = msddaebvp_sol2segs(prob,'',run,lab); % Build 'coll' continuation problem

% T = 2*pi/om2;
dim = 2; N = 5;
uidx_glue = [];
for i=1:2*N+1
fid = sprintf('msbvp.seg%d', i);
prob = coupling_sol2seg(prob, fid, [], run, lab);

fbid = sprintf('msbvp.seg%d.coupling', i);
[data_gi, uidx_gi] = coco_get_func_data(prob, fbid, 'data', 'uidx');

uidx_glue = [uidx_glue; uidx_gi(data_gi.gamma_idx(2))];
end

fbid = coco_get_id('msbvp.seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata1, uidx1] = coco_get_func_data(prob, fbid, 'data', 'uidx');

uidx_glue = [uidx_glue; uidx1(fdata1.p_idx(3)); uidx1(fdata1.T_idx)];
data = [];
data.nsegs = 2*N+1;
prob = coco_add_func(prob,'gluing',@glue_coup,@glue_coup_du,data,'zero',...
                   'uidx', uidx_glue);
               
prob = coco_add_func(prob,'fun1',@fun1,[],'inactive',...
                   'om2','uidx', uidx1(fdata1.T_idx));
               
prob = coco_set(prob,'cont','NPR',1);
% prob = coco_set(prob,'cont','ItMX',500);
prob = coco_set(prob,'cont','atlas','kd','PtMX',40000,...
    'NAdapt',0,'R',.01,'R_max',10,'R_min',1e-4, 'MaxRes', 10,'al_max',15);
% coco(prob,'run_test2',[],1,{'om2','om1','alpha'},{[0.75,1.5]})
coco(prob,'run_surface_new',[],2,{'alpha','om2','om1','Texc'},{[0.1,2],[0.4,2.5],[0,2]})





%% Run to find a fold point
run = 'run_test1'; lab = 5;

prob = coco_prob();

prob = msddaebvp_sol2segs(prob,'',run,lab); % Build 'coll' continuation problem

% T = 2*pi/om2;
dim = 2; N = 5;
uidx_glue = [];
for i=1:2*N+1
fid = sprintf('msbvp.seg%d', i);
prob = coupling_sol2seg(prob, fid, [], run, lab);

fbid = sprintf('msbvp.seg%d.coupling', i);
[data_gi, uidx_gi] = coco_get_func_data(prob, fbid, 'data', 'uidx');

uidx_glue = [uidx_glue; uidx_gi(data_gi.gamma_idx(2))];
end

fbid = coco_get_id('msbvp.seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata1, uidx1] = coco_get_func_data(prob, fbid, 'data', 'uidx');

uidx_glue = [uidx_glue; uidx1(fdata1.p_idx(3)); uidx1(fdata1.T_idx)];
data = [];
data.nsegs = 2*N+1;
prob = coco_add_func(prob,'gluing',@glue_coup,@glue_coup_du,data,'zero',...
                   'uidx', uidx_glue);
% a=1;
prob = coco_set(prob,'cont','MaxRes',10);
prob = coco_set(prob,'cont','ItMX',2000);
prob = coco_set(prob,'cont','NPR',1);
% prob = coco_set(prob,'corr','ItMX',100);
% prob = coco_set(prob,'cont','atlas','kd','PtMX',40000,...
%     'NAdapt',0,'R',.01,'R_max',25,'R_min',1e-4,'al_max',15);

% coco(prob,'run_singlefold',[],1,{'om1','om2','alpha'},{[],[0.4,2.5]})

% coco(prob,'run_test2',[],1,{'om1','alpha','om2'},{[-0.15,1.2],[0.1,5]})

coco(prob,'run_test2',[],1,{'om1','Texc','alpha'},{[],[3.92,12.56]})

% Plotting the fold branch
% bd = coco_bd_read('run_folds_2');
% col = coco_bd_col(bd,{'alpha','om2','om1'});
% hold on
% plot3(col(1,:)',col(2,:)',col(3,:)','b-.','Linewidth',2)



%% Continuation run to find a fold/branch point

run = 'run_test2'; lab = 189;

prob = coco_prob();
%====Constructing the zero problem
%1. Trajectory segments and boundary conditions
prob = msddaebvp_sol2segs(prob,'',run,lab); % Build 'coll' continuation problem

dim = 2; N = 5;
%2. Coupling conditions
uidx_glue = [];
for i=1:2*N+1
fid = sprintf('msbvp.seg%d', i);
prob = coupling_sol2seg(prob, fid, [], run, lab);

fbid = sprintf('msbvp.seg%d.coupling', i);
[data_gi, uidx_gi] = coco_get_func_data(prob, fbid, 'data', 'uidx');

uidx_glue = [uidx_glue; uidx_gi(data_gi.gamma_idx(2))];
end

%3. Gluing conditions
fbid = coco_get_id('msbvp.seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata1, uidx1] = coco_get_func_data(prob, fbid, 'data', 'uidx');
uidx_glue = [uidx_glue; uidx1(fdata1.p_idx(3)); uidx1(fdata1.T_idx)];
data = [];
data.nsegs = 2*N+1;
prob = coco_add_func(prob,'gluing',@glue_coup,@glue_coup_du,data,'zero',...
                   'uidx', uidx_glue);


%====Constructing the Adjoint problem
%1. Adjoint of the trajectory segments and boundary conditions
prob = adjt_msddaebvp_isol2segs(prob,''); % Build 'coll' continuation problem

%2. Adjoints of the coupling conditions
axidx_glue = []; 
for i=1:2*N+1
fid = sprintf('msbvp.seg%d', i);
prob = adjt_coupling_isol2seg(prob, fid, []);

fbid = sprintf('msbvp.seg%d.coupling', i);
[adata_gi, axidx_gi] = coco_get_adjt_data(prob, fbid, 'data', 'axidx');
copti = adata_gi.coupling_opt;

axidx_glue = [axidx_glue; axidx_gi(copti.gamma_adjt_idx(2))];
end

[adata1, axidx1] = coco_get_adjt_data(prob, 'msbvp.seg1.ddaecoll', 'data', 'axidx');
opt1  = adata1.ddaecoll_opt;
axidx_glue = [axidx_glue; axidx1(opt1.p_idx(3)); axidx1(opt1.T_idx)];

prob = coco_add_adjt(prob,'gluing','aidx',axidx_glue);


prob = coco_set(prob,'cont','NPR',1);
prob = coco_set(prob,'cont','MaxRes',10);
prob = coco_set(prob,'cont','hmax',100000);
prob = coco_set(prob,'cont','ItMX',500);
coco(prob,'run_test3',[],1,{'om1','Texc','d.om1','d.alpha'},{[],[4.2,4.24]})
% coco(prob,'run_test4',[],1,{'om1','om2','d.om1','d.alpha'},{[],[1.46,1.47]})
% coco(prob,'run_test3',[],1,{'om1','om2','d.om1','d.alpha'},{[],[1.15,1.27]})

% coco(prob,'run_test_fold2',[],1,{'om1','om2','d.om1','d.alpha'},{[],[0.4,1.6]})


%% Run to drive the continuation parameter to unity

run = 'run_test3'; lab = 10;

prob = coco_prob();
%====Constructing the zero problem
%1. Trajectory segments and boundary conditions
prob = msddaebvp_sol2segs(prob,'',run,lab); % Build 'coll' continuation problem

dim = 2; N = 5;
%2. Coupling conditions
uidx_glue = [];
for i=1:2*N+1
fid = sprintf('msbvp.seg%d', i);
prob = coupling_sol2seg(prob, fid, [], run, lab);

fbid = sprintf('msbvp.seg%d.coupling', i);
[data_gi, uidx_gi] = coco_get_func_data(prob, fbid, 'data', 'uidx');

uidx_glue = [uidx_glue; uidx_gi(data_gi.gamma_idx(2))];
end

%3. Gluing conditions
fbid = coco_get_id('msbvp.seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata1, uidx1] = coco_get_func_data(prob, fbid, 'data', 'uidx');
uidx_glue = [uidx_glue; uidx1(fdata1.p_idx(3)); uidx1(fdata1.T_idx)];
data = [];
data.nsegs = 2*N+1;
prob = coco_add_func(prob,'gluing',@glue_coup,@glue_coup_du,data,'zero',...
                   'uidx', uidx_glue);


%====Constructing the Adjoint problem
%1. Adjoint of the trajectory segments and boundary conditions
prob = adjt_msddaebvp_sol2segs(prob,'',run,lab,'BP'); % Build 'coll' continuation problem

%2. Adjoints of the coupling conditions
axidx_glue = []; 
for i=1:2*N+1
fid = sprintf('msbvp.seg%d', i);
prob = adjt_coupling_sol2seg(prob, fid,[],run,lab,'BP');

fbid = sprintf('msbvp.seg%d.coupling', i);
[adata_gi, axidx_gi] = coco_get_adjt_data(prob, fbid, 'data', 'axidx');
copti = adata_gi.coupling_opt;

axidx_glue = [axidx_glue; axidx_gi(copti.gamma_adjt_idx(2))];
end

[adata1, axidx1] = coco_get_adjt_data(prob, 'msbvp.seg1.ddaecoll', 'data', 'axidx');
opt1  = adata1.ddaecoll_opt;
axidx_glue = [axidx_glue; axidx1(opt1.p_idx(3)); axidx1(opt1.T_idx)];

% branch switch data
chart = coco_read_solution(run, lab, 'chart');
cdata = coco_get_chart_data(chart, 'lsol');

[chart, aidx] = coco_read_adjoint('gluing', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob,'gluing','aidx',axidx_glue,'l0',chart.x,...
                                          'tl0',cdata.v(aidx));
prob = coco_set(prob,'cont','NPR',1);
prob = coco_set(prob,'cont','h_max',1000);
prob = coco_set(prob,'cont','MaxRes',10);
prob = coco_set(prob,'cont','ItMX',500);
% prob = coco_set(prob,'cont','atlas','kd','PtMX',40000,...
%     'NAdapt',0,'R',.01,'R_max',25,'R_min',1e-4,'al_max',10);
% prob = coco_set(prob,'cont','atlas','kd','PtMX',40000,...
%     'NAdapt',0,'R',0.01,'R_min',1e-4,'R_max',25,'al_max',15);

% coco(prob,'run_test5',[],1,{'d.om1','om1','om2','d.alpha'},{[0, 1]})
coco(prob,'run_test4',[],1,{'d.om1','om1','Texc','d.alpha'},{[0, 1]})

%% Driving d.alpha to zero
run = 'run_test4'; lab = 6;
prob = coco_prob();
%====Constructing the zero problem
%1. Trajectory segments and boundary conditions
prob = msddaebvp_sol2segs(prob,'',run,lab); % Build 'coll' continuation problem

dim = 2; N = 5;
%2. Coupling conditions
uidx_glue = [];
for i=1:2*N+1
fid = sprintf('msbvp.seg%d', i);
prob = coupling_sol2seg(prob, fid, [], run, lab);

fbid = sprintf('msbvp.seg%d.coupling', i);
[data_gi, uidx_gi] = coco_get_func_data(prob, fbid, 'data', 'uidx');

uidx_glue = [uidx_glue; uidx_gi(data_gi.gamma_idx(2))];
end

%3. Gluing conditions
fbid = coco_get_id('msbvp.seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata1, uidx1] = coco_get_func_data(prob, fbid, 'data', 'uidx');
uidx_glue = [uidx_glue; uidx1(fdata1.p_idx(3)); uidx1(fdata1.T_idx)];
data = [];
data.nsegs = 2*N+1;
prob = coco_add_func(prob,'gluing',@glue_coup,@glue_coup_du,data,'zero',...
                   'uidx', uidx_glue);


%====Constructing the Adjoint problem
%1. Adjoint of the trajectory segments and boundary conditions
prob = adjt_msddaebvp_sol2segs(prob,'',run,lab); % Build 'coll' continuation problem

%2. Adjoints of the coupling conditions
axidx_glue = []; 
for i=1:2*N+1
fid = sprintf('msbvp.seg%d', i);
prob = adjt_coupling_sol2seg(prob, fid, [], run,lab);

fbid = sprintf('msbvp.seg%d.coupling', i);
[adata_gi, axidx_gi] = coco_get_adjt_data(prob, fbid, 'data', 'axidx');
copti = adata_gi.coupling_opt;

axidx_glue = [axidx_glue; axidx_gi(copti.gamma_adjt_idx(2))];
end

[adata1, axidx1] = coco_get_adjt_data(prob, 'msbvp.seg1.ddaecoll', 'data', 'axidx');
opt1  = adata1.ddaecoll_opt;
axidx_glue = [axidx_glue; axidx1(opt1.p_idx(3)); axidx1(opt1.T_idx)];

% branch switch data
chart = coco_read_solution(run, lab, 'chart');
% cdata = coco_get_chart_data(chart, 'lsol');

[chart, aidx] = coco_read_adjoint('gluing', run, lab, 'chart', 'lidx');
% prob = coco_add_adjt(prob,'gluing','aidx',axidx_glue,'l0',chart.x,...
%                                           'tl0',cdata.v(aidx));
prob = coco_add_adjt(prob,'gluing','aidx',axidx_glue,'l0',chart.x);

prob = coco_set(prob,'cont','h_max',1000);
prob = coco_set(prob,'cont','MaxRes',10);
% coco(prob,'run_test6',[],1,{'d.om2','om1','alpha','om2'},{[],[],[0.4,3],[0.1,2]})
coco(prob,'run_test5',[],1,{'d.alpha','om1','Texc','alpha'},{[],[],[pi,4*pi],[0.1,2.5]})



%% Driving alpha

run = 'run_upperridge_part2'; lab = 42;

prob = coco_prob();
%====Constructing the zero problem
%1. Trajectory segments and boundary conditions
prob = msddaebvp_sol2segs(prob,'',run,lab); % Build 'coll' continuation problem

dim = 2; N = 5;
%2. Coupling conditions
uidx_glue = [];
for i=1:2*N+1
fid = sprintf('msbvp.seg%d', i);
prob = coupling_sol2seg(prob, fid, [], run, lab);

fbid = sprintf('msbvp.seg%d.coupling', i);
[data_gi, uidx_gi] = coco_get_func_data(prob, fbid, 'data', 'uidx');

uidx_glue = [uidx_glue; uidx_gi(data_gi.gamma_idx(2))];
end

%3. Gluing conditions
fbid = coco_get_id('msbvp.seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata1, uidx1] = coco_get_func_data(prob, fbid, 'data', 'uidx');
uidx_glue = [uidx_glue; uidx1(fdata1.p_idx(3)); uidx1(fdata1.T_idx)];
data = [];
data.nsegs = 2*N+1;
prob = coco_add_func(prob,'gluing',@glue_coup,@glue_coup_du,data,'zero',...
                   'uidx', uidx_glue);


%====Constructing the Adjoint problem
%1. Adjoint of the trajectory segments and boundary conditions
prob = adjt_msddaebvp_sol2segs(prob,'',run,lab); % Build 'coll' continuation problem

%2. Adjoints of the coupling conditions
axidx_glue = []; 
for i=1:2*N+1
fid = sprintf('msbvp.seg%d', i);
prob = adjt_coupling_sol2seg(prob, fid, [], run,lab);

fbid = sprintf('msbvp.seg%d.coupling', i);
[adata_gi, axidx_gi] = coco_get_adjt_data(prob, fbid, 'data', 'axidx');
copti = adata_gi.coupling_opt;

axidx_glue = [axidx_glue; axidx_gi(copti.gamma_adjt_idx(2))];
end

[adata1, axidx1] = coco_get_adjt_data(prob, 'msbvp.seg1.ddaecoll', 'data', 'axidx');
opt1  = adata1.ddaecoll_opt;
axidx_glue = [axidx_glue; axidx1(opt1.p_idx(3)); axidx1(opt1.T_idx)];

% branch switch data
chart = coco_read_solution(run, lab, 'chart');
% cdata = coco_get_chart_data(chart, 'lsol');

[chart, aidx] = coco_read_adjoint('gluing', run, lab, 'chart', 'lidx');
% prob = coco_add_adjt(prob,'gluing','aidx',axidx_glue,'l0',chart.x,...
%                                           'tl0',cdata.v(aidx));
prob = coco_add_adjt(prob,'gluing','aidx',axidx_glue,'l0',chart.x);

% prob = coco_set(prob,'cont','ItMX',200);
prob = coco_set(prob,'cont','NPR',1);
prob = coco_set(prob,'cont','atlas','kd','PtMX',40000,...
    'NAdapt',0,'R',.01,'R_max',250,'R_min',1e-4,'al_max',15,'MaxRes',100);

coco(prob,'run_upperridge_part3',[],1,{'d.alpha','om1','om2','alpha'},{[],[],[],[1.5963,2]})


%% Plotting the quasiperiodic surface

% figure: 20a
% xp = 0.3; yp = 2.5; zp = -0.2;
xp = 0.8; yp = 2; zp = -0.2;

%==folds: (1): 1.2,2*pi/4.229,0.202  (2) 1.2,2*pi/4.724,0.181 
%         (3): 1.2,2*pi/6.62,0.305   (4): 0.75,2*pi/3.605,0.773

%====fold branch
bd = coco_bd_read('run_singlefold');
col = coco_bd_col(bd,{'alpha','om2','om1'});
hold on
plot3(col(1,:)',col(2,:)',col(3,:)','b-.','Linewidth',2)

l1 = length(col(2,:)');
hold on
plot3(col(1,:)',col(2,:)',zp*ones(l1),'b-.','Linewidth',1)

hold on
plot3(col(1,:)',yp*ones(l1),col(3,:)','b-.','Linewidth',1)

hold on
plot3(xp*ones(l1),col(2,:)',col(3,:)','b-.','Linewidth',1)

grid('on')
box('on') 
xlim([0,2])

hold on
plot3(1.2,2*pi/4.229,0.202,'b-.','Linewidth',1)


%=======upperridge
bd = coco_bd_read('run_upperridge');
col = coco_bd_col(bd,{'alpha','om2','om1'});
hold on
plot3(col(1,:)',col(2,:)',col(3,:)','k-','Linewidth',2)

l1 = length(col(2,:)');
hold on
plot3(col(1,:)',col(2,:)',zp*ones(l1),'k-','Linewidth',1)

hold on
plot3(col(1,:)',yp*ones(l1),col(3,:)','k-','Linewidth',1)

hold on
plot3(xp*ones(l1),col(2,:)',col(3,:)','k-','Linewidth',1)

%====lowerridge 1:
bd = coco_bd_read('run_lowerridge1');
col = coco_bd_col(bd,{'alpha','om2','om1'});
hold on
plot3(col(1,:)',col(2,:)',col(3,:)','r-','MarkerSize',10)

l1 = length(col(2,:)');
hold on
plot3(col(1,:)',col(2,:)',zp*ones(l1),'r-','MarkerSize',6)

hold on
plot3(col(1,:)',yp*ones(l1),col(3,:)','r-','MarkerSize',6)

hold on
plot3(xp*ones(l1),col(2,:)',col(3,:)','r-','MarkerSize',6)

hold on
% plot3(0.75,2*pi/3.605,0.775,'r.','MarkerSize',30,'MarkerFaceColor','r')
a=plot3(1.2,2*pi/6.62,0.305,'r.','MarkerSize',30,'MarkerFaceColor','r');
uistack(a,'top')

set(gca,'FontSize',12)
xlabel('$\mu_{\alpha}$', 'Interpreter', 'Latex', 'Fontsize',18)
ylabel('$2\pi/\mu_{T}$', 'Interpreter', 'Latex', 'Fontsize',18)
zlabel('$\mu_{\omega}$', 'Interpreter', 'Latex', 'Fontsize',18)
H=gca;
H.LineWidth = 1.5;

xlim([0,1.98])
ylim([0.35,2.5])
zlim([-0.2,1.4])

%==final manifold obtained with varying alpha
bd = coco_bd_read('run_lowerridge2');
col = coco_bd_col(bd,{'alpha','om2','om1'});
hold on
plot3(col(1,:)',col(2,:)',col(3,:)','g-','LineWidth',3)

l1 = length(col(2,:)');
hold on
plot3(col(1,:)',col(2,:)',zp*ones(l1),'g-','Linewidth',1)

hold on
plot3(col(1,:)',yp*ones(l1),col(3,:)','g-','Linewidth',1)

hold on
plot3(xp*ones(l1),col(2,:)',col(3,:)','g-','Linewidth',1)

hold on
a=plot3(1.2,2*pi/4.724,0.181,'k.','MarkerSize',30);
uistack(a,'top')


xlim([0.8,2])
ylim([0.4,2])
zlim([-0.2,1.5])


bd = coco_bd_read('run_lowerridge2');
col = coco_bd_col(bd,{'alpha','om2','om1'});
hold on
plot3(col(1,:)',col(2,:)',col(3,:)','g-','LineWidth',3)


atlas = coco_bd_read('run_surface_new', 'atlas');
plot_atlas_kd(atlas.charts,1,2,3)
grid('on')
box('on')


