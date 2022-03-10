%===condition number and frequency response computations
coco_use_recipes_toolbox ddecoll_vz1
%% Frequency response curve

run='run2';

Amp = []; omega=[]; phase = [];
for i=1:124
 lab = i;
 [sol1,data] = ddaecoll_read_solution('seg1',run,lab);
 w = sol1.p(1);   
 amp = max(abs(sol1.x(:,1)));
 Amp =[Amp;amp];
 
 x = interp1(sol1.t/sol1.t(end),sol1.x,data.tbp);
%  figure(1);
%  plot(sol1.t/sol1.t(end),sol1.x(:,1))
%  hold on
%  plot(data.tbp,x(:,1),'ro')

 po.solution.profile = x';
 po.solution.degree = 4;
 po.solution.tbp = reshape(data.tbp,[5,10]);
 po.solution.lab = i;
 pophas=arrayfun(@(d)coll_phase(d.solution),po);
 
%  omega = [omega;w];
omega = [omega;w];
 phase = [phase;pophas];
end
figure(1);
plot(omega,Amp,'b-','LineWidth',2)

hold on
plot(omega,phase,'r-','LineWidth',2)

set(gca,'Fontsize',12)
xlabel('$$\omega$$','Interpreter','Latex','Fontsize',18)
legend({'$$\mathrm{response\,amplitude}$$',...
       ['$$\mathrm{response\,phase}$$' newline '$$(\mathrm{divided\,by}2\pi)$$']},...
        'Interpreter','Latex','LineWidth',1.5,'Fontsize',14)
legend('boxoff')     
box('on')
grid('on')
H=gca;
H.LineWidth = 1.5;

%% plotting the response
run = 'run4'; lab = 14;
[sol1,data] = ddaecoll_read_solution('seg1',run,lab);
figure(4)
plot(sol1.t/sol1.t(end),sol1.x(:,1))
hold on
plot(sol1.t/sol1.t(end),sol1.x(:,2))
w=sol1.p(2);
sqrt(sol1.x(1,1)^2+(sol1.x(1,2)/w)^2)


%% Condition number computations

run = 'run2';
Cond = []; 
for i=1:124
lab = i;
[sol1,data] = ddaecoll_read_solution('seg1',run,lab);
% w = sol1.p(2); 
% Omega = [Omega;w]; 
 
chart = coco_read_solution(run,lab,'chart');
Jac=full(chart.private.data{1,2}.cond{2});


Jac(104:203,206:207) = Jac(104:203,206:207)/10;

ind = [100,101];
Jac(ind,:)=[];




% cnd=1/min(svds(Jac,1,'smallest','SubspaceDimension',60));
cnd=1/min(svds(Jac,1,'smallest'));
Cond = [Cond; cnd];
end
figure(3);
hold on
plot(phase,log10(5e-3*Cond),'m-','LineWidth',2)

set(gca,'Fontsize',12)
box('on')
grid('on')
xlabel('$$\mathrm{response\,phase\,(divided\,by\,2\pi)}$$','Interpreter','Latex','Fontsize',14)
ylabel('$$\mathrm{log}_{10}\left(\zeta/\mathrm{minimal\,singular\,value} \right)$$','Interpreter','Latex',...
         'Fontsize',14)

H=gca;
H.LineWidth = 1.5;

% legend({'$$\mathrm{fixed}\,a,\,\mathrm{fixed}\,b$$',...
%        '$$\mathrm{fixed}\,b$$','$$\mathrm{fixed}\,\sqrt{a^2+b^2}$$',...
%          '$$\mathrm{released\,a,b}$$'},'Interpreter','Latex','LineWidth',1.5,'Fontsize',14,'NumColumns',2)

legend({'$$\mathrm{fixed}\,a,\,\mathrm{fixed}\,b$$',...
       '$$\mathrm{fixed}\,b$$','$$\mathrm{fixed}\,\sqrt{a^2+b^2}$$',...
         '$$\mathrm{full\,3d\,manifold}$$'},'Interpreter','Latex','LineWidth',1.5,'Fontsize',14,'NumColumns',2)
     
legend('boxoff')

% xlim([0.05,0.5])


%% ===== 


Jac1 = A29;
Jac1 = Jac1([100,101],:);
cnd1 = 1/min(svds(Jac1,1,'smallest'))


Jac2 = A93;

for j=206:207
for i=104:203
    Jac2(i,j) = Jac2(i,j)/10;
end
end
Jac2([100,101],:)=[];
cnd2 = 1/min(svds(Jac2,1,'smallest'))












