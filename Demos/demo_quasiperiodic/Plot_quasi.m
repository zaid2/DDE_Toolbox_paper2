%Plotting the quasiperiodic case

bd = coco_bd_read('run_multiplefolds');

for i=1
  figure(i+1); 
  clf; hold on; grid on
  
%    [sol, data] = msddebvp_read_solution('', 'torus_run1', 7);
  [sol, data] = msddaebvp_read_solution('', 'run_multiplefolds', 55);
  N  = data.nsegs;
  M  = size(sol{1}.x,1);
  x0 = zeros(N+1,2);
  x1 = zeros(N+1,2);
  XX = zeros(M,N+1);
  YY = XX;
  ZZ = XX;
  for j=1:N+1
    n       = mod(j-1,N)+1;
    XX(:,j) = sol{n}.x(1:M,1);
    ZZ(:,j) = sol{n}.x(1:M,2);
    YY(:,j) = sol{n}.t(1:M);
    x0(j,:) = sol{n}.x(1,:);
    x1(j,:) = sol{n}.x(M,:);
     sol{j}.T = sol{1}.t(end);
  end
  
  surf(XX,YY/sol{1}.t(end),ZZ, 'FaceColor', 0.9*[1 1 1], 'FaceAlpha', 0.7, ...
    'MeshStyle', 'column', 'LineStyle', '-', 'EdgeColor', 0.6*[1 1 1], ...
    'LineWidth', 2);
  plot3(x0(:,1), zeros(N+1,1), x0(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', 'black', 'Marker', '.', 'MarkerSize', 12);
  plot3(x1(:,1), ones(N+1,1)*sol{1}.t(end)/sol{1}.t(end), x1(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', 'black', 'Marker', '.', 'MarkerSize', 12);
  
  hold off; view([50 15])
  box('on')
  xlabel('$x_{1}$','Interpreter','Latex','Fontsize',18)
  ylabel('$\tau$','Interpreter','Latex','Fontsize',18)
  zlabel('$x_{2}$','Interpreter','Latex','Fontsize',18)
  
end


%=====plot adjoints
for i=1:11
fid = sprintf('msbvp.seg%d', i);

[sol,data] = coupling_read_adjoint(fid,[],'run_test4',7);

% plot3(sol.tbp,sol.l0_coupling(1,:),sol.l0_coupling(2,:))
plot(sol.l0_coupling(1,:),sol.l0_coupling(2,:))
hold on
% fbid = sprintf('msbvp.seg%d.coupling', i);
% [adata_gi, axidx_gi] = coco_get_adjt_data(prob, fbid, 'data', 'axidx');
% copti = adata_gi.coupling_opt;
% 
% axidx_glue = [axidx_glue; axidx_gi(copti.Delta_adjt_idx);...
%              axidx_gi(copti.gamma_adjt_idx)];
end

