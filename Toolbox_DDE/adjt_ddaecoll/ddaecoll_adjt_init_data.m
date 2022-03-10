function data = ddaecoll_adjt_init_data(prob, src_data) %#ok<INUSL>
%COLL_ADJT_INIT_DATA   Initialize data structure for a 'coll' adjoint problem.
%
% See also ODE_INIT_DATA, COCO_SAVE_DATA, COCO_FUNC_DATA.

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: ode_init_data.m 2902 2015-10-09 18:06:32Z hdankowicz $

data.ddaecoll_seg = src_data;

seg  = src_data.ddaecoll;
NCOL = seg.NCOL;
NTST = seg.NTST;
dim  = src_data.dim;
ydim = src_data.ydim;
pdim = src_data.pdim;
xbp_shp = src_data.xbp_shp;
tbp_idx = src_data.tbp_idx;


bpdim  = dim*(NCOL+1);      % Number of basepoint values per interval
xbpdim = dim*(NCOL+1)*NTST; % Number of basepoint values
ybpdim = ydim*(NCOL+1)*NTST;    % Number of collocation conditions for y
cndim  = dim*NCOL;          % Number of collocation node values per interval
xcnnum = NCOL*NTST;         % Number of collocation nodes
xcndim = dim*NCOL*NTST;     % Number of collocation conditions             
cntnum = NTST-1;            % Number of internal boundaries
cntdim = dim*(NTST-1);      % Number of continuity conditions

addim  = xcndim+cntdim+2*dim+NTST*(NCOL+1)*ydim+2+pdim;

opt.adim    = [xbpdim, addim];
opt.xcn_idx = (1:xcndim)';
opt.x0_idx  = xcndim + cntdim + (1:dim)';
opt.x1_idx  = opt.x0_idx(end) + (1:dim)';
opt.ybp_idx = opt.x1_idx(end) + (1:ybpdim)';
opt.T0_idx  = opt.ybp_idx(end) + 1;
opt.T_idx   = opt.T0_idx + 1;
opt.p_idx   = opt.T_idx + (1:pdim)';

opt.fbp_idx = 1:xbpdim;

opt.id0 = [eye(dim); zeros(xbpdim-dim,dim)];
opt.id1 = [zeros(xbpdim-dim,dim); eye(dim)];

% rows = 1:(NTST*NCOL+1)*dim;
% cols = reshape(1:xbpdim, xbp_shp);
% cols = cols(:, tbp_idx);

rows = 1:xbpdim;
rows = reshape(rows,[dim,NTST*(NCOL+1)]);
opt.dybprows = repmat(rows,[ydim,1]);
opt.dybpcols = repmat(1:NTST*(NCOL+1)*ydim,[dim,1]);

% opt.Wbp = sparse(rows,cols,ones((NTST*NCOL+1)*dim,1));

tcn         = src_data.taucn;
opt.dTtcn   = repmat(tcn', [dim 1]);

data.ddaecoll_opt = opt;

data = init_data(data);
end



function data = init_data(data)
%INIT_DATA   Initialize data for COLL adjoint problem.

seg  = data.ddaecoll_seg;
opt  = data.ddaecoll_opt;

NCOL = seg.ddaecoll.NCOL;
NTST = seg.ddaecoll.NTST;
dim = seg.dim;
pdim = seg.pdim;
ydim = seg.ydim;

cndim  = NCOL*dim;
bpdim  = dim*(NCOL+1);
xbpdim = NTST*(NCOL+1)*dim;
xcnnum = NTST*NCOL;
xcndim = NTST*NCOL*dim;
ybpdim = NTST*(NCOL+1)*ydim;
ycndim = NTST*NCOL*ydim;
cntdim = (NTST-1)*dim;

addim  = xcndim+cntdim+2*dim+ybpdim+2+pdim;


opt.dJrows = xcndim;
opt.dJcols = addim*(xbpdim+ybpdim+2+pdim);

opt.dpdTtcn = repmat(permute(seg.taucn, [2 3 1]), [dim, pdim]);
opt.dTtcn   = repmat(seg.taucn', [dim 1]);
opt.dxdTtcn = repmat(permute(seg.taucn, [2 3 1]), [dim, dim]);
opt.dydTtcn = repmat(permute(seg.taucn, [2 3 1]), [dim, ydim]);

tbp = seg.tbp;
opt.dydTtbp = repmat(permute(tbp, [2 3 1]), [dim,ydim]);



% for Jacobian of adjoint with respect to delta_x
rows = repmat(1:dim, [1 dim]);
rows = repmat(rows(:), [1 dim]) + xcndim*repmat(0:dim-1, [dim^2 1]);
rows = repmat(rows(:), [1 xcnnum]) + dim*repmat(0:xcnnum-1, [dim^3 1]);
opt.dxdxrows1 = rows;
cols = kron(1:xcndim, ones(1,dim));
opt.dxdxcols1 = repmat(reshape(cols, [dim^2 xcnnum]), [dim 1]);

dxdxrows = repmat(reshape(1:xcndim, [cndim NTST]), [dim*bpdim, 1]);
cols = 1 + dim*(0:NCOL-1);
cols = repmat(cols(:), [1 dim]) + repmat(0:dim-1, [NCOL 1]);
cols = repmat(cols(:), [1 bpdim]) + addim*repmat(0:bpdim-1, [cndim 1]);
cols = repmat(cols(:), [1 NTST]) + ...
  (cndim+addim*bpdim)*repmat(0:NTST-1, [cndim*bpdim 1]);
dxdxcols = kron(cols, ones(dim,1));

idx = 1:cndim;
idx = repmat(idx(:), [1 dim]) + xcndim*repmat(0:dim-1, [cndim 1]);
idx = repmat(idx(:), [1 bpdim]) + ...
  dim*xcndim*repmat(0:bpdim-1, [dim*cndim 1]);
idx = repmat(idx(:), [1 NTST]) + ...
  dim*(NCOL+cndim*xbpdim)*repmat(0:NTST-1, [dim*cndim*bpdim 1]);
opt.dxdxidx = idx(:);

rows = reshape(1:dim*dim*NCOL*NTST,[dim*dim xcnnum]);
opt.dxdyrows1 = repmat(rows, [ydim, 1]);
opt.dxdycols1 = repmat(1:ycndim, [dim*dim, 1]);

rows = reshape(1:xcndim, [dim,xcnnum]);
dxdyrows = repmat(rows, [ydim, xbpdim]);
cols = xbpdim*addim+repmat(1:xcndim,[dim,1]);
dxdycols = repmat(cols(:),[1,ybpdim])+ repmat(0:addim:(ybpdim-1)*addim,[dim*xcndim, 1]);

dxdT0rows = repmat(reshape(1:xcndim, [dim xcnnum]), [dim 1]);
dxdT0cols = (xbpdim+ybpdim)*addim + repmat(1:xcndim, [dim 1]);
dxdTrows  = repmat(reshape(1:xcndim, [dim xcnnum]), [dim 1]);
dxdTcols  = (xbpdim+ybpdim+1)*addim + repmat(1:xcndim, [dim 1]);

if pdim~=0
dxdprows = repmat(reshape(1:xcndim, [dim xcnnum]), [dim*pdim 1]);
cols = (xbpdim+ybpdim+2)*addim + repmat(1:dim, [dim 1]);
cols = repmat(cols(:), [1 pdim]) + addim*repmat(0:pdim-1, [dim^2 1]);
cols = repmat(cols(:), [1 xcnnum]) + ...
  dim*repmat(0:xcnnum-1, [dim^2*pdim 1]);
dxdpcols = cols;
else
dxdprows = []; dxdpcols = [];    
end

opt.dxrows = [dxdxrows(:); dxdyrows(:); dxdT0rows(:); dxdTrows(:); dxdprows(:)];
opt.dxcols = [dxdxcols(:); dxdycols(:); dxdT0cols(:); dxdTcols(:); dxdpcols(:)];

% for Jacobian of adjoint with respect to delta_y

rows = reshape(1:xbpdim, [dim,NTST*(NCOL+1)]);
dydxrows = repmat(rows,[ydim*dim,1]);
cols = reshape(1:ybpdim,[ydim,NTST*(NCOL+1)]);
cols = repelem(cols,1,dim)+repmat(0:addim:(xbpdim-1)*addim,[ydim,1]);
dydxcols = xcndim+cntdim+2*dim+repmat(cols(:)',[dim,1]);


rows = reshape(1:xbpdim, [dim,NTST*(NCOL+1)]);
dydyrows  = repmat(rows,[ydim*ydim,1]);
cols = reshape(1:ybpdim,[ydim,NTST*(NCOL+1)]);
cols = repelem(cols,1,ydim)+repmat(0:addim:(ybpdim-1)*addim,[ydim,1]);
dydycols = xbpdim*addim+xcndim+cntdim+2*dim+repmat(cols(:)',[dim,1]);



rows = reshape(1:xbpdim, [dim,NTST*(NCOL+1)]);
cols = 1:NTST*(NCOL+1)*ydim;
dydT0rows = repmat(rows, [ydim,1]);
dydT0cols = (xbpdim+ybpdim)*addim+xcndim+cntdim+2*dim+repmat(cols, [dim,1]);
dydTrows  = repmat(rows, [ydim,1]);
dydTcols  = (xbpdim+ybpdim+1)*addim+xcndim+cntdim+2*dim+repmat(cols, [dim,1]);


rows = reshape(1:xbpdim, [dim,NTST*(NCOL+1)]);
dydprows  = repmat(rows,[ydim*pdim,1]);

cols = reshape(1:NTST*(NCOL+1)*ydim, [ydim, NTST*(NCOL+1)]);
cols = repelem(cols,1,pdim)+repmat(0:addim:(pdim-1)*addim, [ydim, NTST*(NCOL+1)]);
dydpcols = (xbpdim+ybpdim+2)*addim+xcndim+cntdim+2*dim+repmat(cols(:)',[dim,1]);

opt.dyrows = [dydxrows(:); dydyrows(:); dydT0rows(:); dydTrows(:); dydprows(:)];
opt.dycols = [dydxcols(:); dydycols(:); dydT0cols(:); dydTcols(:); dydpcols(:)];

% for Jacobian of adjoint with respect to delta_T0 & delta_T
dT0dxrows = repmat(reshape(1:xcndim, [cndim NTST]), [bpdim 1]);
dT0dxcols = addim-pdim-1+addim*repmat(0:xbpdim-1, [cndim 1]);

idx = 1:cndim;
idx = repmat(idx(:), [1 bpdim]) + xcndim*repmat(0:bpdim-1, [cndim 1]);
idx = repmat(idx(:), [1 NTST]) + ...
  (cndim+bpdim*xcndim)*repmat(0:NTST-1, [cndim*bpdim 1]);
opt.dT0dxidx = idx(:);
opt.dTdxidx  = opt.dT0dxidx;

dT0dyrows = repmat((1:xcndim)',[1,ybpdim]);
cols = xbpdim*addim+xcndim+cntdim+2*dim+ybpdim+ones(xcndim,ybpdim);
dT0dycols = cols + repmat(0:addim:(ybpdim-1)*addim, [xcndim,1]);

dT0dT0rows = 1:xcndim;
dT0dTrows  = 1:xcndim;
dT0dT0cols = (xcndim+cntdim+2*dim+ybpdim+1+(xbpdim+ybpdim)*addim)*ones(xcndim,1);
dT0dTcols  = (xcndim+cntdim+2*dim+ybpdim+1+(xbpdim+ybpdim+1)*addim)*ones(xcndim,1);

if pdim~=0
dT0dprows = repmat(reshape(1:xcndim, [dim xcnnum]), [pdim 1]);
cols = xcndim+cntdim+2*dim+ybpdim+(xbpdim+ybpdim+2)*addim + ones(dim,1);
cols = repmat(cols(:), [1 pdim]) + addim*repmat(0:pdim-1, [dim 1]);
dT0dpcols = repmat(cols(:), [1 xcnnum]);
else
dT0dprows=[]; dT0dpcols = [];    
end
dT0rows = [dT0dxrows(:); dT0dyrows(:); dT0dT0rows(:); dT0dTrows(:); dT0dprows(:)];
dT0cols = [dT0dxcols(:); dT0dycols(:); dT0dT0cols(:); dT0dTcols(:); dT0dpcols(:)];

% for Jacobian of adjoint with respect to delta_p

if pdim~=0
rows = reshape(1:xcndim*pdim, [dim xcnnum pdim]);
rows = permute(repmat(rows, [dim 1 1]), [1 3 2]);
opt.dpdxrows1 = rows(:,:);
cols = kron(1:xcndim, ones(1,dim));
opt.dpdxcols1 = repmat(reshape(cols, [dim^2 xcnnum]), [pdim 1]);

dpdxrows2 = repmat(reshape(1:xcndim, [cndim NTST]), [pdim*bpdim 1]);
cols = repmat(1:pdim, [cndim, 1]);
dpdxcols2 = repmat(cols(:), [1, xbpdim]) + ...
  xcndim+cntdim+2*dim+ybpdim+2+addim*repmat(0:xbpdim-1, [cndim*pdim 1]);

idx = 1:cndim;
idx = repmat(idx(:), [1 pdim*bpdim]) + ...
  xcndim*repmat(0:pdim*bpdim-1, [cndim 1]);
idx = repmat(idx(:), [1 NTST]) + ...
  (cndim+pdim*bpdim*xcndim)*repmat(0:NTST-1, [cndim*pdim*bpdim 1]);
opt.dpdxidx = idx(:);


rows = reshape(1:xcndim,[dim,xcnnum]);
rows = repmat(rows, [ydim*pdim, 1]);
opt.dpdyrows1 = rows(:);


cols = (reshape(1:ycndim*pdim, [ycndim,pdim]))';
cols = repelem(cols, dim,1);
opt.dpdycols1 = cols(:);

rows = repmat((1:xcndim)',[1,ybpdim*pdim]);
dpdyrows = rows(:);
cols = xbpdim*addim +xcndim+cntdim+2*dim+ybpdim+2+repmat(1:addim:ybpdim*addim,[xcndim,1]);
cols = repmat(cols, [1,pdim]) + repelem(0:pdim-1,xcndim,ybpdim);
dpdycols = cols(:);

dpdT0rows = seg.dprows;
dpdT0cols = seg.dpcols + xcndim+cntdim+2*dim+ybpdim+2+(xbpdim+ybpdim)*addim;
dpdTrows  = seg.dprows;
dpdTcols  = seg.dpcols + xcndim+cntdim+2*dim+ybpdim+2+(xbpdim+ybpdim+1)*addim;

dpdprows = repmat(reshape(1:xcndim, [dim, xcnnum]), [pdim^2, 1]);
cols = repmat(1:pdim, [dim 1]);
cols = repmat(cols(:), [1 pdim]) + xcndim+cntdim+2*dim+ybpdim+2 + ...
  (xbpdim+ybpdim+2)*addim+addim*repmat(0:pdim-1, [dim*pdim 1]);
cols = repmat(cols(:), [1 xcnnum]);
dpdpcols = cols;

dprows = [dpdxrows2(:); dpdyrows(:); dpdT0rows(:); dpdTrows(:); dpdprows(:)];
dpcols = [dpdxcols2(:); dpdycols(:); dpdT0cols(:); dpdTcols(:); dpdpcols(:)];

else
dprows = [];
dpcols = [];
opt.dpdxidx = [];
end


opt.dT0Tprows = [dT0rows; dT0rows; dprows];
opt.dT0Tpcols = [dT0cols; 1+dT0cols; dpcols];

data.ddaecoll_opt  = opt;


end


