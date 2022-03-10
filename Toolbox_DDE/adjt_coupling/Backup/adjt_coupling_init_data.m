function data = adjt_coupling_init_data(prob, src_data) %#ok<INUSL>

seg  = src_data.ddaecoll_seg;
data.coupling_seg = src_data;
coupling = data.coupling_seg; % coupling function data

NTST = seg.ddaecoll.NTST; % Number of mesh intervals
NCOL = seg.ddaecoll.NCOL; % Degree of polynomial interpolants
dim  = seg.dim;           % State-space dimension
pdim = seg.pdim;          % Number of problem parameters
ydim = seg.ydim;
xbpdim = seg.xbpdim;
yiid = coupling.yiid;

% yccndim  = dim*NCOL*NTST;
xcndim  = seg.xcndim;     % Number of collocation conditions
cntdim = dim*(NTST-1);

xi = coupling.xi;
Ji = coupling.Ji;
T_uniq  = coupling.T_uniq;
Ci = coupling.Ci;

T_idx = coupling.T_idx;
Delta_idx = coupling.Delta_idx;
gamma_idx = coupling.gamma_idx;

addim   = xcndim*length(Ji)+(NTST*NCOL+1)*dim+length(T_idx)+pdim+...
                                length(Delta_idx)+length(gamma_idx)+cntdim;


opt.addim  = addim;
opt.adim   = [xbpdim addim];

eqns_adjt = xcndim+cntdim+2*dim+(NTST*NCOL+1)*ydim+2+pdim;

xcn_eqidx = (1:xcndim)' + (Ji(:)-1)'*eqns_adjt;
temp = repmat(yiid(:),[1,(NTST*NCOL+1)])+repmat(0:ydim:ydim*(NTST*NCOL), [dim,1]);
yibp_eqidx = xcndim+cntdim+2*dim+temp(:)'+eqns_adjt*(xi-1);
T_eqidx   = xcndim+cntdim+2*dim+(NTST*NCOL+1)*ydim+2+eqns_adjt*(T_uniq(:)-1);
p_eqidx   = xcndim+cntdim+2*dim+(NTST*NCOL+1)*ydim+2+(1:pdim)';

opt.guidx_adjt        = [xcn_eqidx(:); yibp_eqidx(:); T_eqidx; p_eqidx];

opt.xcn_adjt_idx   = (1:xcndim*length(Ji))';
opt.yibp_adjt_idx  = opt.xcn_adjt_idx(end)+(1:(NTST*NCOL+1)*dim)';
opt.T_adjt_idx     = opt.yibp_adjt_idx(end)+(1:length(T_eqidx))';
opt.p_adjt_idx     = opt.T_adjt_idx(end)+(1:pdim)';
opt.Delta_adjt_idx = length(opt.guidx_adjt) + (1:length(Delta_idx))';
opt.gamma_adjt_idx = length(opt.guidx_adjt)+length(Delta_idx)+...
                                                    (1:length(gamma_idx))';
% if ~isempty(data.coupling_seg.gamma_idx)
% 
% end

wts = seg.wts;
opt.wts = wts;

opt.tcn       = seg.M;
                                            
XBPDIM   = length(coupling.Xibp_idx);
YBPDIM   = length(coupling.yibp_idx);
TDIM     = length(coupling.T_idx);
Deltadim = length(coupling.Delta_idx);
Gammadim = length(coupling.gamma_idx);

opt.dJrows = xbpdim;
opt.dJcols = addim*(XBPDIM+YBPDIM+TDIM+Deltadim+pdim+Gammadim); %====cols = addim*Nunknowns

data.coupling_opt = opt;


end



