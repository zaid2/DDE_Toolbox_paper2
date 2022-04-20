function data = adjt_coupling_init_data(prob, src_data) %#ok<INUSL>
% Copyright (C) Zaid Ahsan, Mingwu Li
seg  = src_data.ddaecoll_seg;
data.coupling_seg = src_data;
coupling_seg = data.coupling_seg; % coupling function data

NTST = seg.ddaecoll.NTST; % Number of mesh intervals
NCOL = seg.ddaecoll.NCOL; % Degree of polynomial interpolants
dim  = seg.dim;           % State-space dimension
pseg_dim = seg.pdim;          % Number of problem parameters
ydim = seg.ydim;
xbpdim = seg.xbpdim;
xcndim  = seg.xcndim;     % Number of collocation conditions
ybpdim = seg.ybpdim;
cntdim = dim*(NTST-1);

yiid  = coupling_seg.yiid;
xi    = coupling_seg.xi;
Ji    = coupling_seg.Ji;
T_uniq  = coupling_seg.T_uniq;
Nj = coupling_seg.Nj;     % Total xj involved in the copling conditions
Ci = coupling_seg.Ci;

T_idx = coupling_seg.T_idx;
Delta_idx = coupling_seg.Delta_idx;
gamma_idx = coupling_seg.gamma_idx;

pdim = coupling_seg.pdim;
dim_xiek = coupling_seg.dim_xiek;
dim_xibk = coupling_seg.dim_xibk;

addim   = xcndim*Nj+NTST*(NCOL+1)*dim+length(T_idx)+Nj*2*dim+pdim+...
                                         length(Delta_idx)+length(gamma_idx); 
                                     %==eqns = coll. nodes+base points(y)+T+pseg+x0+x1+Delta+gamma+pnew+gamma_mesh

opt.addim  = addim;
opt.adim   = [xbpdim+(Ci+1)+(Ci-1)*dim+dim_xiek+dim_xibk addim]; % Lagrange mult. has length = xbpdim+(Ci-1)*dim


%% Equations from the coupling conditions
eqns_adjt = xcndim+cntdim+2*dim+NTST*(NCOL+1)*ydim+2+pseg_dim;

if Nj~=0
xcn_eqidx = (1:xcndim)' + (Ji(:)-1)'*eqns_adjt;
else
xcn_eqidx = 0;    % Only history seg. involved
end

temp = repmat(yiid(:),[1,NTST*(NCOL+1)])+repmat(0:ydim:(ybpdim-1), [dim,1]);
yibp_eqidx = xcndim+cntdim+2*dim+temp(:)'+eqns_adjt*(xi-1);
T_eqidx   = xcndim+cntdim+2*dim+NTST*(NCOL+1)*ydim+2+eqns_adjt*(T_uniq(:)-1);
pseg_eqidx   = xcndim+cntdim+2*dim+NTST*(NCOL+1)*ydim+2+(1:pseg_dim)';

%% Equations from the interval boundary conditions
x0_eqidx = data.coupling_seg.ddaecoll_opt.x0_idx;
x1_eqidx = data.coupling_seg.ddaecoll_opt.x1_idx;

if Nj~=0
   X0_eqidx = x0_eqidx(:)+repmat((Ji(:)-1)'*eqns_adjt,[dim,1]); 
   X1_eqidx = x1_eqidx(:)+repmat((Ji(:)-1)'*eqns_adjt,[dim,1]); 
   
   opt.guidx_adjt        = [xcn_eqidx(:); yibp_eqidx(:); T_eqidx; X0_eqidx(:); X1_eqidx(:); pseg_eqidx;];
else
   opt.guidx_adjt        = [yibp_eqidx(:); T_eqidx; pseg_eqidx]; % No history segment involved
   X0_eqidx = [];
end



if Nj~=0
opt.xcn_adjt_idx   = (1:xcndim*Nj)';
opt.yibp_adjt_idx  = opt.xcn_adjt_idx(end)+(1:NTST*(NCOL+1)*dim)';
else
opt.yibp_adjt_idx  = (1:NTST*(NCOL+1)*dim)';    
end
opt.T_adjt_idx     = opt.yibp_adjt_idx(end)+(1:length(T_eqidx))';
if ~isempty(X0_eqidx)                                                
opt.X0_adjt_idx = opt.T_adjt_idx(end) + (1:length(X0_eqidx(:)))';
opt.X1_adjt_idx = opt.T_adjt_idx(end) + length(X0_eqidx(:)) + (1:length(X1_eqidx(:)))';
opt.p_adjt_idx     = opt.X1_adjt_idx(end)+(1:pdim)';
else
opt.X0_adjt_idx = [];    
opt.X1_adjt_idx = [];  
opt.p_adjt_idx     = opt.T_adjt_idx(end)+(1:pdim)';
end

opt.Delta_adjt_idx = opt.p_adjt_idx(end) + (1:length(Delta_idx))';
opt.gamma_adjt_idx = opt.Delta_adjt_idx(end)+(1:length(gamma_idx))';

pnew_dim = coupling_seg.pnew_dim;
if pnew_dim~=0
opt.pnew_idx = opt.T_adjt_idx + pseg_dim + (1:pnew_dim)';
else
opt.pnew_idx = [];    
end

wts = seg.wts;
opt.wts = wts;

opt.tcn       = seg.taucn;

if Nj~=0
XBPDIM   = length(coupling_seg.Xibp_idx);
else
XBPDIM = 0;    
end
YBPDIM   = length(coupling_seg.yibp_idx);
TDIM     = length(coupling_seg.T_idx);
Deltadim = length(coupling_seg.Delta_idx);
Gammadim = length(coupling_seg.gamma_idx);

opt.dJrows = xbpdim+(Ci+1)+(Ci-1)*dim+dim_xiek+dim_xibk;
opt.dJcols = addim*(XBPDIM+YBPDIM+TDIM+pdim+Deltadim+Gammadim); %====cols = addim*Nunknowns

data.coupling_opt = opt;


end



