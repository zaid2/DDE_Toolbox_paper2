%% Duffing oscillator: Frequency Response Curve

zeta  = 5*10^(-3);
mu    = 1;
a     = 7.5*10^(-3);
b     = 0;
gamma = -0.01;
alpha = 1;
Omega = 0.5;
T     = 2*pi/Omega;

[t0,x0,y0]=init_guess_duffing(alpha,Omega,gamma,mu,zeta,a,b);

prob = coco_prob();
prob = coco_set(prob, 'ddaecoll', 'NTST', 10);


prob = ddaecoll_isol2seg(prob,'seg1',@duffing,@duffing_dfdt,@duffing_dfdx,...
    @duffing_dfdy,@duffing_dfdp,...
    t0,x0,y0,{'Omega','a','b'},...
    [Omega,a,b]);

% Adding the boundary conditions
fbid = coco_get_id('seg1', 'ddaecoll'); % Create toolbox instance identifier
[fdata, uidx] = coco_get_func_data(prob, fbid, 'data', 'uidx');

prob = coco_add_func(prob,'bc',@bc, fdata, 'zero','uidx',...
    [uidx(fdata.T0_idx);...
    uidx(fdata.T_idx); uidx(fdata.p_idx(1))]);

T = 2*pi/Omega;

% Coupling conditions
g1 = {1, 1, -(1-alpha/T),[0,alpha/T]};
g2 = {1, 1, (alpha/T),[alpha/T,1]};
xi = 1;
prob = coupling_isol2seg(prob,'seg1',[],xi,1:2,{g1,g2});                                 

[data_g1, uidx_g1] = coco_get_func_data(prob, 'seg1.coupling', 'data', 'uidx');

% Adding the gluing conditions
data=[];
prob = coco_add_func(prob,'gluing',@glue_coup,data,'zero','uidx',...
    [uidx_g1(data_g1.gamma_idx(2));...
     uidx(fdata.T_idx)]);
                                
% Solving the system 
prob = coco_set(prob, 'cont', 'ItMX', 300);
prob = coco_set(prob,'cont','NPR',1);
prob = coco_set(prob,'lsol','cond',1);
bd = coco(prob, 'run1', [], 1, {'Omega','a','b'},{[0.5,2]});



%% Plotting the frequency response curve

run='run1';
labs = coco_bd_labs(bd);
Amp = []; omega=[];
for i=1:length(labs)
 lab = i;
 [sol1,data] = ddaecoll_read_solution('seg1',run,lab);
 w = sol1.p(1);   
 amp = max(abs(sol1.x(:,1)));
 Amp =[Amp;amp];
 omega = [omega;w];
end
 
figure(1);
plot(omega,Amp,'b-','LineWidth',2)
set(gca,'Fontsize',12)
xlabel('\omega','Fontsize',14)
box('on')
grid('on')
H=gca;
H.LineWidth = 1.5;
 
 