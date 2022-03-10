%% ===PRC of Mackey-Glass Equation

% Run1: To get family of periodic orbits

a=2; b=10; 
alpha=0.4708; 
T=2*pi/sqrt(15); 

omega = 2*pi/T;

prob = coco_prob();

t0 = (0:T/100:T)';
x0 = 0.01*cos(omega*t0-pi/2)+1;
y0 = 0.01*cos(omega*(t0-alpha)-pi/2)+1;

prob = coco_set(prob, 'ddaecoll', 'NTST', 10);
prob = ddaecoll_isol2seg(prob,'seg1',@mackey_f,@mackey_dfdt,@mackey_dfdx,...
    @mackey_dfdy,@mackey_dfdp,@mackey_dfdtdt,@mackey_dfdtdx,...
    @mackey_dfdtdy,@mackey_dfdtdp,@mackey_dfdxdx,@mackey_dfdxdy,...
    @mackey_dfdxdp,@mackey_dfdydy,@mackey_dfdydp,@mackey_dfdpdp,...
    t0,x0,y0,{'a','b','alpha'},...
    [a,b,alpha]);

% Adding the boundary conditions
fbid = coco_get_id('seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata, uidx] = coco_get_func_data(prob, fbid, 'data', 'uidx');

prob = coco_add_func(prob,'bc',@bc,@bc_du, fdata, 'zero','uidx',...
    [uidx(fdata.x0_idx); uidx(fdata.T0_idx)]);

% Coupling conditions
g1 = {1, 1, -(1-alpha/T),[0,alpha/T]};
% g1 = {j_indices, Delta, [gamma_ikb,gamma_ike]}
g2 = {1, 1, (alpha/T),[alpha/T,1]};
xi = 1;
prob = coupling_isol2seg(prob,'seg1',[],xi,1,{g1,g2});                                 

[data_g1, uidx_g1] = coco_get_func_data(prob, 'seg1.coupling', 'data', 'uidx');

% Adding the gluing conditions
data=[];
prob = coco_add_func(prob,'gluing',@glue_coup,@glue_coup_du,data,'zero','uidx',...
    [uidx_g1(data_g1.gamma_idx(2)); uidx(fdata.T_idx); uidx(fdata.p_idx(3))]);
                                 
prob = coco_add_pars(prob, 'obj', fdata.T_idx, 'T');

prob = coco_set(prob,'corr','ItMX',20);

% Solving the system 
prob = coco_set(prob, 'cont', 'h_max', 10);
prob = coco_set(prob, 'cont', 'ItMX', 100);
prob = coco_set(prob, 'cont', 'NPR', 4);
% bd = coco(prob, 'run1', [], 1, {'T' 'alpha'},{[],[0.4708,1.1]});
bd1 = coco(prob, 'run1', [], 1, {'T' 'alpha'},{[],[0.4708,1.1]});

%% Run2: Phase response curve of specific orbit

prob = coco_prob();
lab=6; run='run1';

prob = coco_set(prob, 'ddaecoll', 'NTST', 10);
prob = ddaecoll_sol2seg(prob, 'seg1', run, lab); % Reconstruct 'coll' continuation problem
prob = adjt_ddaecoll_isol2seg(prob, 'seg1');

% Adding the boundary conditions
fbid = coco_get_id('seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata, uidx] = coco_get_func_data(prob, fbid, 'data', 'uidx');
prob = coco_add_func(prob,'bc',@bc,@bc_du, fdata, 'zero','uidx',...
       [uidx(fdata.x0_idx); uidx(fdata.T0_idx)]);

[data, axidx] = coco_get_adjt_data(prob, 'seg1.ddaecoll', 'data', 'axidx');
opt  = data.ddaecoll_opt;
prob = coco_add_adjt(prob, 'bc', 'aidx', ...
   axidx([opt.x0_idx; opt.T0_idx]));

% Coupling conditions
prob = coupling_sol2seg(prob, 'seg1', [], run, lab);
[data_g1, uidx_g1] = coco_get_func_data(prob, 'seg1.coupling', 'data', 'uidx');

% Adjoint of Coupling Conditions
prob = adjt_coupling_isol2seg(prob, 'seg1', []);

% Adding the gluing conditions
prob = coco_add_func(prob,'gluing',@glue_coup,@glue_coup_du,data,'zero','uidx',...
    [uidx_g1(data_g1.gamma_idx(2)); uidx(fdata.T_idx); uidx(fdata.p_idx(3))]);

[cdata1, cxidx1] = coco_get_adjt_data(prob, 'seg1.coupling', 'data', 'axidx');
copt1 = cdata1.coupling_opt; 
axidx_glue = [cxidx1(copt1.gamma_adjt_idx(2)); axidx(opt.T_idx); axidx(opt.p_idx(3))]; 
prob = coco_add_adjt(prob, 'gluing', 'aidx', axidx_glue); 

prob = coco_add_pars(prob, 'obj', fdata.T_idx, 'T');
prob = coco_add_adjt(prob, 'obj', 'd.T', 'aidx', axidx(opt.T_idx));

axidx = coco_get_adjt_data( prob, 'obj', 'afidx');
fcn  = @(f) @(p,d,u,l,v) deal(d, f(u,l,v));
prob = coco_add_comp(prob, 'fun1', fcn(@(u,l,v) l-1), [], 'zero', ...
   'lidx', axidx);

prob = coco_set(prob, 'corr', 'ItMX', 30);
prob = coco_set(prob, 'cont', 'ItMX', 210);
prob = coco_set(prob,'cont','h_max',1);
% prob = coco_set(prob,'cont','NPR',1);
% bd = coco(prob, 'run2', [], 1, {'b','T','d.a','d.b','d.T'}, {[10,13]});
bd = coco(prob, 'run2', [], 1, {'b' 'T','d.a','d.b','d.T','d.alpha'}, {[9,14]});


%% plotting the solution

%==figure(1): Sample orbits of limit cycles
a=2; b=10;
run='run1';
labs = coco_bd_labs(bd1);
figure(1);
for i=1:length(labs)
lab = i;
[sol1,data]=ddaecoll_read_solution('seg1',run,lab);
chart = coco_read_solution(run,lab,'chart');
x = sol1.x; y = sol1.y;
xdot = a*y./(1+y.^b)-x;
% plot(sol1.t/sol1.t(end),sol1.x-1,'Linewidth',2)
hold on
if i==1
    plot(x,xdot,'LineWidth',3)
else
plot(x,xdot,'LineWidth',2)
end
end
grid on
box('on')
set(gca,'FontSize',12)
xlabel('$x$', 'Interpreter', 'Latex', 'Fontsize',18)
ylabel('$\dot{x}$', 'Interpreter', 'Latex', 'Fontsize',18)
H = gca;
H.LineWidth = 1.5;


%==figure:2 and figure:3
run='run2'; 
lab = 4;
[sol1,data]=ddaecoll_read_solution('seg1',run,lab);
sol2=ddaecoll_read_adjoint('seg1',run,lab);
figure(2);
plot(sol1.t/sol1.t(end),sol1.x,'LineWidth',2)
xlabel('$\tau$', 'Interpreter', 'Latex', 'Fontsize',18)
ylabel('$x$', 'Interpreter', 'Latex', 'Fontsize',18)
grid 'on'
box('on')
H = gca;
H.LineWidth = 1.5;

figure(3);
plot(sol2.tbp,sol2.l0,'LineWidth',2)
set(gca,'FontSize',12)
xlabel('$\tau$', 'Interpreter', 'Latex', 'Fontsize',18)
ylabel('$\lambda_{\mathrm{DE}}$', 'Interpreter', 'Latex', 'Fontsize',18)
grid 'on'
box('on')
H = gca;
H.LineWidth = 1.5;


% fig:3 and fig:4
% periodic orbits
% run='plot_data'; 

% [sol1,data]=ddaecoll_read_solution('seg1',run,1);
% hold on
% plot(sol1.t/sol1.t(end),sol1.x)

% x = sol1.t/sol1.t(end);
% y = 9:0.5:14;
% [X,Y] = meshgrid(x,y);
% Z1 = zeros(size(X,1),size(X,2));
% Z2 = zeros(size(X,1),size(X,2));
% 
% for i=1:length(y)
% lab = i;    
% [sol1,data]=ddaecoll_read_solution('seg1',run,lab);
% sol2=ddaecoll_read_adjoint('seg1',run,lab);
% 
% Z1(i,:) = sol1.x;
% Z2(i,:) = sol2.l0;
% box('on')
% grid('on')
% 
% end
% 
% figure(3);
% surf(X',Y',Z1', 'FaceColor', 0.9*[1 1 1], 'FaceAlpha', 1,'MeshStyle', 'column','LineStyle', '-','EdgeColor', 0.6*[1 1 1], ...
%     'LineWidth', 0.5);
% set(gca,'FontSize',12)
% xlabel('$\tau$', 'Interpreter', 'Latex', 'Fontsize',18)
% ylabel('$b$', 'Interpreter', 'Latex', 'Fontsize',18)
% ylabel('$x$', 'Interpreter', 'Latex', 'Fontsize',18)
% grid 'on'
% box('on')
% H = gca;
% H.LineWidth = 1.5;
% zlim([0.63,1.25])
% ylim([9,14])
% 
% 
% figure(4);
% surf(X',Y',Z2', 'FaceColor', 0.9*[1 1 1], 'FaceAlpha', 1,'MeshStyle', 'column','LineStyle', '-','EdgeColor', 0.6*[1 1 1], ...
%     'LineWidth', 0.5);
% set(gca,'FontSize',12)
% xlabel('$\tau$', 'Interpreter', 'Latex', 'Fontsize',18)
% ylabel('$b$', 'Interpreter', 'Latex', 'Fontsize',18)
% zlabel('$\lambda_{\mathrm{DE}}$', 'Interpreter', 'Latex', 'Fontsize',18)
% grid on
% box('on')
% H = gca;
% H.LineWidth = 1.5;
% zlim([-2.5,2.5])
% ylim([9,14])










