%====function plot function

%=====plotting alpha vs om2 vs om1

bd = coco_bd_read('run_test6');
om1   = coco_bd_col(bd,'om1');
om2   = coco_bd_col(bd,'om2');
alpha = coco_bd_col(bd,'alpha');
hold on
% plot3(alpha,om2,om1,'Linewidth',2)
plot3(alpha,om2,om1,'r.')
xlim([0.6,1.4])


% labs = coco_bd_labs(bd);
for i=1
  figure(i+1); 
  clf; hold on; grid on
  
%    [sol, data] = msddebvp_read_solution('', 'torus_run1', 7);
  [sol, data] = msddaebvp_read_solution('', 'run6', 3);
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
  
  surf(XX,YY,ZZ, 'FaceColor', 0.9*[1 1 1], 'FaceAlpha', 0.7, ...
    'MeshStyle', 'column', 'LineStyle', '-', 'EdgeColor', 0.6*[1 1 1], ...
    'LineWidth', 0.5);
  plot3(x0(:,1), zeros(N+1,1), x0(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', 'black', 'Marker', '.', 'MarkerSize', 12);
  plot3(x1(:,1), ones(N+1,1)*sol{1}.t(end), x1(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', 'black', 'Marker', '.', 'MarkerSize', 12);
  
  hold off; view([50 15])
  box('on')
  xlabel('$x_{1}$','Interpreter','Latex','Fontsize',15)
  ylabel('$t$','Interpreter','Latex','Fontsize',15)
  zlabel('$x_{2}$','Interpreter','Latex','Fontsize',15)
  
end


%=====plot adjoints
for i=1:11
fid = sprintf('msbvp.seg%d', i);

[sol,data] = coupling_read_adjoint(fid,'run6',3);

plot3(sol.tbp,sol.l0(1,:),sol.l0(2,:))
hold on
% fbid = sprintf('msbvp.seg%d.coupling', i);
% [adata_gi, axidx_gi] = coco_get_adjt_data(prob, fbid, 'data', 'axidx');
% copti = adata_gi.coupling_opt;
% 
% axidx_glue = [axidx_glue; axidx_gi(copti.Delta_adjt_idx);...
%              axidx_gi(copti.gamma_adjt_idx)];
end


