%% Run to find the zero solution from the Hopf point

% Copyright (C) Zaid Ahsan

alpha=0.8255; p1=0.5; p2 = -1.257;
T = 2*pi/0.868;
t0 = (0:0.01:T)';
omega = 0.868;
 r=0.01;
x0 = [1+r*omega^-1*sin(omega*t0-pi/2) r*cos(omega*t0-pi/2)];
y0 = [1+r*omega^-1*sin(omega*(t0-alpha)-pi/2) r*cos(omega*(t0-alpha)-pi/2)];

prob = coco_prob();
prob = coco_set(prob, 'ddaecoll', 'NTST', 100);
prob = ddaecoll_isol2seg(prob,'seg1',@connecting_f,@connecting_dfdt,@connecting_dfdx,...
    @connecting_dfdy,@connecting_dfdp, t0,x0,y0,{'p1','p2','alpha'},...
    [p1,p2,alpha]);

% Adding the boundary conditions
fbid = coco_get_id('seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata, uidx] = coco_get_func_data(prob, fbid, 'data', 'uidx');

% Coupling conditions
g1 = {1, 1, -(1-alpha/T),[0,alpha/T]};
% g1 = {j_indices, Delta, [gamma_ikb,gamma_ike]}
g2 = {1, 1, (alpha/T),[alpha/T,1]};
xi = 1;
prob = coupling_isol2seg(prob,'seg1',[],xi,1:2,{g1,g2});    
[data_g1, uidx_g1] = coco_get_func_data(prob, 'seg1.coupling', 'data', 'uidx');

% Adding the gluing conditions
data=[];
prob = coco_add_func(prob,'gluing',@glue_coup,data,'zero','uidx',...
    [uidx_g1(data_g1.gamma_idx(2)); uidx(fdata.T_idx); uidx(fdata.p_idx(3))]);

prob = coco_add_func(prob,'bc',@bc_per, fdata, 'zero','uidx',...
    [uidx(fdata.T0_idx);uidx(fdata.xbp_idx(2:2:end))],'u0',0.5);

prob = coco_add_pars(prob,'Tcr',2012,'Tcr');


prob = coco_add_pars(prob,'T',uidx(fdata.T_idx),'T');

prob = coco_add_event(prob,'Tpo','T',20);

prob = coco_set(prob, 'cont', 'ItMX', [0,200]);
bd = coco(prob, 'run_per', [], 1,{'p2' 'T' 'Tcr'}); 

%% Connecting orbit: Run to compute p2 and eps

run='run_per'; lab = 11;
[sol1,data]=ddaecoll_read_solution('seg1',run,lab);

chart = coco_read_solution(run, lab, 'chart');

%==scaling the time
scale = 1;
T  = sol1.t(end);
t0 = [sol1.t(1:end-1,1); T*(scale-1)+sol1.t(end,1)];
    
x0 = sol1.x; 
y0 = sol1.y; p1 = sol1.p(1); p2 = sol1.p(2); alpha = sol1.p(3);
% t0 = sol1.t;

prob = coco_prob();
prob = coco_set(prob, 'ddaecoll', 'NTST', 100);
prob = ddaecoll_isol2seg(prob,'seg1',@connecting_f,@connecting_dfdt,@connecting_dfdx,...
    @connecting_dfdy,@connecting_dfdp, t0,x0,y0,{'p1','p2','alpha'},...
    [p1,p2,alpha]);


% Coupling conditions
g1 = {0, 0, alpha,[0,alpha/T],@hist,@hist_dt,@hist_dp,[],[],[]};
% g1 = {j_indices, Delta, [gamma_ikb,gamma_ike]}
g2 = {1, 1, (alpha/T),[alpha/T,1]};
xi = 1;
prob = coupling_isol2seg(prob,'seg1',[],xi,1:2,{g1,g2},{'lambda' 'v1' 'v2' 'eps'},[0.597;0.8585;0.5126;1e-3]);    
[data_g1, uidx_g1] = coco_get_func_data(prob, 'seg1.coupling', 'data', 'uidx');

% Adding the boundary conditions
fbid = coco_get_id('seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata, uidx, u0] = coco_get_func_data(prob, fbid, 'data', 'uidx', 'u0');

prob = coco_add_func(prob,'bc',@bc_conn, fdata, 'zero','uidx',...
    [uidx(fdata.T0_idx);uidx(fdata.x1_idx);uidx(fdata.p_idx(2));...
     uidx_g1(data_g1.p_idx(4:6));uidx(fdata.xbp_idx(2:2:end,1))],'u0',chart.x(end));
 
 % Adding the gluing conditions
data=[];
prob = coco_add_func(prob,'gluing',@glue_coup,data,'zero','uidx',...
    [uidx_g1(data_g1.gamma_idx(2)); uidx(fdata.T_idx);uidx(fdata.p_idx(3))]);

prob = coco_add_pars(prob,'T',uidx(fdata.T_idx),'T');
idx = coco_get_func_data(prob,'bc','uidx');
prob = coco_add_pars(prob,'Tcr',idx(end),'Tcr');

prob = coco_set(prob,'cont','ItMX',[0,50]);
bd = coco(prob, 'run1_homo', [], 0,{'p2' 'eps' 'lambda' 'v1' 'v2'},{[]}); 

%% 2-dimensional continuation of homoclinic orbits

% run='run1_homo'; lab = 1;
% [sol1,data]=ddaecoll_read_solution('seg1',run,lab);
% chart = coco_read_solution(run, lab, 'chart');
% 
% prob = coco_prob();
% prob = coco_set(prob, 'ddaecoll', 'NTST', 100);
% prob = ddaecoll_sol2seg(prob,'seg1',run,lab);
% 
% % % Coupling conditions
% prob = coupling_sol2seg(prob,'seg1',[],run,lab);
% [data_g1, uidx_g1] = coco_get_func_data(prob, 'seg1.coupling', 'data', 'uidx');
% 
% % Adding the boundary conditions
% fbid = coco_get_id('seg1', 'ddaecoll'); % Create toolbox instance identifier
% [fdata, uidx] = coco_get_func_data(prob, fbid, 'data', 'uidx');
% 
% prob = coco_add_func(prob,'bc',@bc_conn, fdata, 'zero','uidx',...
%     [uidx(fdata.T0_idx);uidx(fdata.x1_idx);uidx(fdata.p_idx(2));...
%      uidx_g1(data_g1.p_idx(4:6));uidx(fdata.xbp_idx(2:2:end,1))],'u0',chart.p(end));
% 
%  % Adding the gluing conditions
% data=[];
% prob = coco_add_func(prob,'gluing',@glue_coup,data,'zero','uidx',...
%     [uidx_g1(data_g1.gamma_idx(2)); uidx(fdata.T_idx);uidx(fdata.p_idx(3))]);
% 
% prob = coco_add_pars(prob,'T',uidx(fdata.T_idx),'T');
% idx = coco_get_func_data(prob,'bc','uidx');
% prob = coco_add_pars(prob,'Tcr',idx(end),'Tcr');
% 
% prob = coco_set(prob,'cont','atlas','kd','PtMX',40000,...
%     'NAdapt',0,'R',5,'R_max',5,'R_min',1e-4, 'MaxRes', 5,'al_max',15);
% 
% bd = coco(prob, 'run2_homo', [], 2,{'p1' 'p2' 'alpha' 'T' 'lambda' 'v1' 'v2'},{[0.2,1.4],[],[0.5,1.5]});



%% figures

% figure:1
run = 'run_per';
figure(1);
for i=1:9
lab = i;
[sol1,data]=ddaecoll_read_solution('seg1',run,lab);
figure(1)
hold on
plot(sol1.x(:,1),sol1.x(:,2),'Linewidth',2)
grid on
end
box on
set(gca,'FontSize',12)
title('Family of periodic orbits')
xlabel('$x_{1}$', 'Interpreter', 'Latex', 'Fontsize',18)
ylabel('$x_{2}$', 'Interpreter', 'Latex', 'Fontsize',18)
H = gca;
H.LineWidth = 1.5;

%==figure:2
bd = coco_bd_read('run_per');
col = coco_bd_col(bd,{'p2','T'});
figure(2)
plot(col(1,:)',col(2,:)','b.','Linewidth',2)
grid('on')
box('on')
set(gca,'FontSize',12)
xlabel('$p_{2}$', 'Interpreter', 'Latex', 'Fontsize',18)
ylabel('$T$', 'Interpreter', 'Latex', 'Fontsize',18)
H = gca;
H.LineWidth = 1.5;
xlim([-1.28,-1.035])

%===constant line
y = (7:0.5:62)';
x = -1.0782*ones(length(y),1);
hold on
plot(x,y,'k--')
ylim([0,70])

hold on
% plot(-1.2566082230283877,-5,'rx')
plot(-1.0782,-5,'rx')


% figure 3
run='run1_homo'; lab = 1;
[sol1,data]=ddaecoll_read_solution('seg1',run,lab);
figure(3)
plot(sol1.t/sol1.t(end),sol1.x,'Linewidth',2)
grid on
set(gca,'FontSize',12)
title('Connecting orbit')
xlabel('$\tau$', 'Interpreter', 'Latex', 'Fontsize',18)
ylabel('$x$', 'Interpreter', 'Latex', 'Fontsize',18)
H = gca;
H.LineWidth = 1.5;

% figure:4
% atlas = coco_bd_read('run2_homo', 'atlas');
% plot_atlas_kd(atlas.charts,2,1,3)
% grid('on')
% box('on')
% 
% set(gca,'FontSize',12)
% xlabel('$p_{2}$', 'Interpreter', 'Latex', 'Fontsize',18)
% ylabel('$p_{1}$', 'Interpreter', 'Latex', 'Fontsize',18)
% zlabel('$\alpha$', 'Interpreter', 'Latex', 'Fontsize',18)
% 
% xlim([-2,-0.5])
% ylim([0.46,1])
% zlim([0.5,1.5])
% 
% H = gca;
% H.LineWidth = 1.5;

