function [prob, data] = ddaecoll_construct_adjt(prob, tbid, data, sol)
% Copyright (C) Zaid Ahsan, Mingwu Li

% [data, sol] = init_data(data, sol);
opt  = data.ddaecoll_opt;
seg  = data.ddaecoll_seg;
if ~isempty(sol.l0)
  sol.l0  = interp1(sol.tbp', sol.l0', seg.tbp', 'pchip')';
  sol.l0  = sol.l0(:);
  if ~isempty(sol.tl0)
    sol.tl0 = interp1(sol.tbp', sol.tl0', seg.tbp', 'pchip')';
    sol.tl0 = sol.tl0(:);
  end
end
prob = coco_add_adjt(prob, tbid, @adj, @adj_DU, data, 'l0', sol.l0, ...
  'tl0', sol.tl0, 'adim', opt.adim);
% prob = coco_add_adjt(prob, tbid, @adj, data, 'l0', sol.l0, ...
%   'tl0', sol.tl0, 'adim', opt.adim);
if ~isempty(seg.pnames)
  pfid   = coco_get_id(tbid, 'pars');
  dnames = coco_get_id('d', seg.pnames);
  axidx  = coco_get_adjt_data(prob, tbid, 'axidx');
  prob   = coco_add_adjt(prob, pfid, dnames, 'aidx', axidx(opt.p_idx), ...
    'l0', sol.pars_l0, 'tl0', sol.pars_tl0);
end

end

function [data, J] = adj(prob, data, u) %#ok<INUSL>

seg  = data.ddaecoll_seg;
opt  = data.ddaecoll_opt;
NTST = seg.ddaecoll.NTST;

x  = u(seg.xbp_idx);
y  = u(seg.ybp_idx);
T0 = u(seg.T0_idx);
T  = u(seg.T_idx);

if seg.pdim~=0
p  = u(seg.p_idx);
pcn = repmat(p, seg.p_rep);
pbp  = repmat(p, seg.pbp_rep);
else
pcn = [];
pbp = [];
end

xcn = reshape(seg.W*x, seg.x_shp); % Values at collocation nodes
ycn = reshape(seg.Wy*y, seg.y_shp); 
tcn = seg.taucn*T+T0;

xbp  = reshape(x,seg.xbp_shp);  % Evlaution of variables at N(m+1) base points
ybp  = reshape(y,seg.ybp_shp);
tbp  = seg.tbp*T+T0;

if ~isempty(pcn)
fcn   = seg.fhan(tcn', xcn, ycn, pcn);
fdtcn = seg.dfdthan(tcn', xcn, ycn, pcn);
fdxcn = seg.dfdxhan(tcn', xcn, ycn, pcn);
fdpcn = seg.dfdphan(tcn', xcn, ycn, pcn);

fdybp = seg.dfdyhan(tbp', xbp, ybp, pbp);

else
fcn   = seg.fhan(tcn', xcn, ycn);
fdtcn = seg.dfdthan(tcn', xcn, ycn);
fdxcn = seg.dfdxhan(tcn', xcn, ycn);  

fdybp = seg.dfdyhan(tbp', xbp, ybp);
end

 
% adjoint with respect to delta_x
dxode = sparse(seg.dxrows, seg.dxcols, fdxcn(:));
J = -2*NTST*seg.Wp'-T*seg.W'*dxode;

% adjoint with respect to delta_x^(j+1)(-1), delta_x^(1)(-1), and delta_x^(N)(1)
J = [ J, seg.Q', -opt.id0, opt.id1 ];

% adjoint with respect to delta_ybp
dyode = sparse(opt.dybprows, opt.dybpcols, fdybp(:));
J = [ J, -T*dyode ];


% adjoint with respect to T0 and T
dT0ode = T*fdtcn;
dTode  = fcn + T*fdtcn.*opt.dTtcn;
J = [ J, -(0.5/NTST)*seg.W'*seg.wts2*[dT0ode(:) dTode(:)] ];

% adjoint with respect to p
if seg.pdim~=0
dpode = fdpcn;
dpode = sparse(seg.dprows, seg.dpcols, dpode(:));
J = [ J, -(0.5*T/NTST)*seg.W'*seg.wts2*dpode ];
else
J = [J, []];    
end


end

function [data, dJ] = adj_DU(prob, data, u) %#ok<INUSL>

% [data, dJ1] = coco_ezDFDX('f(o,d,x)', prob, data, @adj, u);

% data = init_data(prob, data);

seg  = data.ddaecoll_seg;
opt  = data.ddaecoll_opt;

NTST = seg.ddaecoll.NTST;
NCOL = seg.ddaecoll.NCOL;
dim  = seg.dim;
ydim = seg.ydim;

x  = u(seg.xbp_idx);
y  = u(seg.ybp_idx);
T0 = u(seg.T0_idx);
T  = u(seg.T_idx);

if seg.pdim~=0
p  = u(seg.p_idx);
pcn = repmat(p, seg.p_rep);
pbp  = repmat(p, seg.pbp_rep);
else
pcn = [];
pbp = [];
end

xcn = reshape(seg.W*x, seg.x_shp);
ycn = reshape(seg.Wy*y, seg.y_shp);
tcn = T0+T*seg.taucn';

xbp  = reshape(x,seg.xbp_shp);  % Evlaution of variables at N(m+1) base points
ybp  = reshape(y,seg.ybp_shp);
tbp  = T0 + seg.tbp'*T;

if ~isempty(pcn)
fcn   = seg.fhan(tcn, xcn, ycn, pcn);
fdtcn = seg.dfdthan(tcn, xcn, ycn, pcn);
fdxcn = seg.dfdxhan(tcn, xcn, ycn, pcn);
fdycn = seg.dfdyhan(tcn, xcn, ycn, pcn);
fdpcn = seg.dfdphan(tcn, xcn, ycn, pcn);
fdybp = seg.dfdyhan(tbp, xbp, ybp, pbp);
else
fcn   = seg.fhan(tcn, xcn, ycn);
fdtcn = seg.dfdthan(tcn, xcn, ycn);
fdxcn = seg.dfdxhan(tcn, xcn, ycn);
fdycn = seg.dfdyhan(tcn, xcn, ycn);
fdybp = seg.dfdyhan(tbp, xbp, ybp);
end

fdxdxcn =   het_DFDXDX(seg, tcn, xcn, ycn, pcn);
fdxdycn =   het_DFDXDY(seg, tcn, xcn, ycn, pcn);
if seg.pdim~=0
fdxdpcn =   het_DFDXDP(seg, tcn, xcn, ycn, pcn);
fdydpcn =   het_DFDYDP(seg, tcn, xcn, ycn, pcn);
fdpdpcn =   het_DFDPDP(seg, tcn, xcn, ycn, pcn);    
fdtdpcn =   het_DFDTDP(seg, tcn, xcn, ycn, pcn);
else
fdxdpcn = []; fdydpcn = []; fdpdpcn = []; fdtdpcn = [];    
end
fdtdtcn =   het_DFDTDT(seg, tcn, xcn, ycn, pcn);
fdtdxcn =   het_DFDTDX(seg, tcn, xcn, ycn, pcn);
fdtdycn =   het_DFDTDY(seg, tcn, xcn, ycn, pcn);

fdxdybp =  het_DFDXDY(seg, tbp, xbp, ybp, pbp);
fdydybp =  het_DFDYDY(seg, tbp, xbp, ybp, pbp); 
fdtdybp =  het_DFDTDY(seg, tbp, xbp, ybp, pbp);
if seg.pdim~=0
fdydpbp =   het_DFDYDP(seg, tbp, xbp, ybp, pbp);
else
fdydpbp = [];    
end

dJrows = opt.dJrows;
dJcols = opt.dJcols;

% Jacobians of adjoint with respect to delta_x
dxdxode  = T*fdxdxcn;
dxdxode  = sparse(opt.dxdxrows1, opt.dxdxcols1, dxdxode(:))*seg.W;
dxdyode  = T*fdxdycn;
dxdyode  = sparse(opt.dxdyrows1, opt.dxdycols1, dxdyode(:))*seg.Wy;
dxdT0ode = T*fdtdxcn;
dxdTode  = fdxcn+T*fdtdxcn.*opt.dxdTtcn;
if seg.pdim~=0
dxdpode  = T*fdxdpcn;
else
dxdpode = [];    
end
dxvals = [dxdxode(opt.dxdxidx); dxdyode(:); dxdT0ode(:); dxdTode(:); dxdpode(:)];

% Jacobain of adjoint with respect to y
fdydxbp = T*permute(fdxdybp,[1 3 2 4]); 
dydxode = fdydxbp;
dydyode = T*fdydybp;
dydT0ode = T*fdtdybp;
dydTode  = fdybp+T*fdtdybp.*opt.dydTtbp;
if seg.pdim~=0
dydpode  = T*fdydpbp;
else
dydpode = [];    
end

dyvals = [dydxode(:); dydyode(:); dydT0ode(:); dydTode(:); dydpode(:)];

% Jacobians of adjoint with respect to T0, T, and p
dT0dxode  = T*fdtdxcn;
dT0dxode  = sparse(seg.dxrows, seg.dxcols, dT0dxode(:))*seg.W;
dT0dyode  = T*fdtdycn;
dT0dyode  = sparse(seg.dyrows, seg.dycols, dT0dyode(:))*seg.Wy;
dT0dT0ode = T*fdtdtcn;
dT0dTode  = fdtcn+T*fdtdtcn.*opt.dTtcn;
if seg.pdim~=0
dT0dpode  = T*fdtdpcn;
else
dT0dpode  = [];    
end
dTdxode   = fdxcn+T*fdtdxcn.*opt.dxdTtcn;
dTdxode   = sparse(seg.dxrows, seg.dxcols, dTdxode(:))*seg.W;
dTdyode   = fdycn+T*fdtdycn.*opt.dydTtcn;
dTdyode   = sparse(seg.dyrows, seg.dycols, dTdyode(:))*seg.Wy;
dTdT0ode  = fdtcn+T*fdtdtcn.*opt.dTtcn;
dTdTode   = (2*fdtcn+T*fdtdtcn.*opt.dTtcn).*opt.dTtcn;
if seg.pdim~=0
dTdpode   = fdpcn+T*fdtdpcn.*opt.dpdTtcn;
else
dTdpode =[];    
end

if seg.pdim~=0
dpdxode  =  T*fdxdpcn;
dpdxode  = sparse(opt.dpdxrows1, opt.dpdxcols1, dpdxode(:))*seg.W;
dpdyode  = T*permute(fdydpcn,[1 3 2 4]);
dpdyode  = sparse(opt.dpdyrows1, opt.dpdycols1, dpdyode(:))*(kron(eye(seg.pdim),seg.Wy));
dpdT0ode = T*fdtdpcn;
dpdTode  = fdpcn+T*fdtdpcn.*opt.dpdTtcn;
dpdpode  = T*fdpdpcn;
else
dpdxode = []; dpdyode=[]; dpdT0ode=[]; dpdTode=[]; dpdpode=[];    
end

dT0vals = [dT0dxode(opt.dT0dxidx); dT0dyode(:); dT0dT0ode(:); dT0dTode(:); dT0dpode(:)];
dTvals  = [dTdxode(opt.dTdxidx); dTdyode(:); dTdT0ode(:); dTdTode(:); dTdpode(:)];
dpvals  = [dpdxode(opt.dpdxidx); dpdyode(:); dpdT0ode(:); dpdTode(:); dpdpode(:)];

dT0Tpvals = [dT0vals; dTvals; dpvals];

dJ = sparse(opt.dT0Tprows, opt.dT0Tpcols, dT0Tpvals, dJrows, dJcols);
dJ =  -(0.5/seg.ddaecoll.NTST)*seg.W'*seg.wts2*dJ;
dJ = dJ - sparse(opt.dyrows, opt.dycols, dyvals, NTST*(NCOL+1)*dim, dJcols);
dJ = dJ - seg.W'*sparse(opt.dxrows, opt.dxcols, dxvals, dJrows, dJcols);

% dJ_incorr = reshape(full(dJ), [opt.adim,length(u)]);
% [data, dJ_corr] = coco_ezDFDX('f(o,d,x)', prob, data, @adj, u);
% check = dJ_corr-dJ_incorr;
% max(abs(check(:)))
end



function Jxx = het_DFDXDX(data, t, x, y, p)
%DFDXDX   Vectorized evaluation of non-autonomous d2F/dx2.

xdim = data.dim;
ydim = data.ydim;

if ~isempty(p)
    if isempty(data.dfdxdxhan)
        f  = @(x,p) data.dfdxhan(p(ydim+1,:), x, p(1:ydim,:), p(ydim+2:end,:));
        ytp = [ y; t ; p ];
        Jxx = coco_ezDFDX('f(x,p)', f, x, ytp);
    else
        [m, n] = size(x);
        Jxx = zeros(m,m,m,n);
        dfdxdxhan = data.dfdxdxhan;
        for i=1:n
            Jxx(:,:,:,i) = dfdxdxhan(t(i), x(:,i), y(:,i), p(:,i));
        end
    end
    
else
    
    if isempty(data.dfdxdxhan)
        f  = @(x,p) data.dfdxhan(p(ydim+1,:), x, p(1:ydim,:));
        yt = [ y; t ];
        Jxx = coco_ezDFDX('f(x,p)', f, x, yt);
    else
        [m, n] = size(x);
        Jxx = zeros(m,m,m,n);
        dfdxdxhan = data.dfdxdxhan;
        for i=1:n
            Jxx(:,:,:,i) = dfdxdxhan(t(i), x(:,i), y(:,i));
        end
    end
  
end

end

function Jxy = het_DFDXDY(data, t, x, y, p)
%DFDXDX   Vectorized evaluation of non-autonomous d2F/dx2.

xdim = data.dim;
ydim = data.ydim;

if ~isempty(p)
    if isempty(data.dfdxdyhan)
        f  = @(x,p) data.dfdxhan(p(xdim+1,:), p(1:xdim,:), x, p(xdim+2:end,:));
        xtp = [ x ; t; p ];
        Jxy = coco_ezDFDX('f(x,p)', f, y, xtp);
    else
        [m, n] = size(x);
        o = size(y,1);
        Jxy = zeros(m,m,o,n);
        dfdxdyhan = data.dfdxdyhan;
        for i=1:n
            Jxy(:,:,:,i) = dfdxdyhan(t(i), x(:,i), y(:,i), p(:,i));
        end
    end
    
else
    
    if isempty(data.dfdxdyhan)
        f  = @(x,p) data.dfdxhan(p(xdim+1,:), p(1:xdim,:), x);
        xt = [ x ; t];
        Jxy = coco_ezDFDX('f(x,p)', f, y, xt);
    else
        [m, n] = size(x);
        o = size(y,1);
        Jxy = zeros(m,m,o,n);
        dfdxdyhan = data.dfdxdyhan;
        for i=1:n
            Jxy(:,:,:,i) = dfdxdyhan(t(i), x(:,i), y(:,i));
        end
    end
    
end
end


function Jxp = het_DFDXDP(data, t, x, y, p)
%DFDXDP   Vectorized evaluation of non-autonomous d2F/dxdp.

xdim = data.dim;
ydim = data.ydim;

  if isempty(data.dfdxdphan)
    f  = @(x,p) data.dfdxhan(x(end,:), x(1:xdim,:), x(xdim+1:xdim+ydim,:), p);
    xyt = [x; y; t ];
    Jxp = coco_ezDFDP('f(x,p)', f, xyt, p);
  else
    [m, n] = size(x);
    o  = size(p,1);
    Jxp = zeros(m,m,o,n);
    dfdxdphan = data.dfdxdphan;
    for i=1:n
      Jxp(:,:,:,i) = dfdxdphan(t(i), x(:,i), y(:,i), p(:,i));
    end
  end

end

function Jyy = het_DFDYDY(data, t, x, y, p)
%DFDXDX   Vectorized evaluation of non-autonomous d2F/dx2.

xdim = data.dim;
ydim = data.ydim;

if ~isempty(p)
    if isempty(data.dfdydyhan)
        f  = @(x,p) data.dfdyhan(p(xdim+1,:), p(1:xdim,:), x, p(xdim+2:end,:));
        xtp = [ x; t ; p ];
        Jyy = coco_ezDFDX('f(x,p)', f, y, xtp);
    else
        [m, n] = size(y);
        o = size(x,1);
        Jyy = zeros(o,m,m,n);
        dfdydyhan = data.dfdydyhan;
        for i=1:n
            Jyy(:,:,:,i) = dfdydyhan(t(i), x(:,i), y(:,i), p(:,i));
        end
    end

else
    if isempty(data.dfdydyhan)
        f  = @(x,p) data.dfdyhan(p(xdim+1,:), p(1:xdim,:), x);
        xt = [ x; t ];
        Jyy = coco_ezDFDX('f(x,p)', f, x, xt);
    else
        [m, n] = size(y);
        o = size(x,1);
        Jyy = zeros(o,m,m,n);
        dfdydyhan = data.dfdydyhan;
        for i=1:n
            Jyy(:,:,:,i) = dfdydyhan(t(i), x(:,i), y(:,i));
        end
    end
  
end
end


function Jyp = het_DFDYDP(data, t, x, y, p)
%DFDXDP   Vectorized evaluation of non-autonomous d2F/dxdp.

xdim = data.dim;
ydim = data.ydim;

  if isempty(data.dfdydphan)
    f  = @(x,p) data.dfdyhan(x(end,:), x(1:xdim,:), x(xdim+1:xdim+ydim,:), p);
    xyt = [x; y; t ];
    Jyp = coco_ezDFDP('f(x,p)', f, xyt, p);
  else
    [m, n] = size(y);
    
    o  = size(p,1);
    Jyp = zeros(size(x,1),m,o,n);
    dfdydphan = data.dfdydphan;
    for i=1:n
      Jyp(:,:,:,i) = dfdydphan(t(i), x(:,i), y(:,i), p(:,i));
    end
  end

end





function Jpp = het_DFDPDP(data, t, x, y, p)
%DFDXDP   Vectorized evaluation of non-autonomous d2F/dxdp.

xdim = data.dim;
ydim = data.ydim;
  if isempty(data.dfdpdphan)
    f  = @(x,p) data.dfdphan(x(end,:), x(1:xdim,:), x(xdim+1:xdim+ydim,:), p);
    xyt = [ x; y; t ];
    Jpp = coco_ezDFDP('f(x,p)', f, xyt, p);
  else
    [m, n] = size(x);
    o  = size(p,1);
    Jpp = zeros(m,o,o,n);
    dfdpdphan = data.dfdpdphan;
    for i=1:n
      Jpp(:,:,:,i) = dfdpdphan(t(i), x(:,i), y(:,i), p(:,i));
    end
  end
end

function Jtx = het_DFDTDX(data, t, x, y, p)
%DFDTDX   Vectorized evaluation of non-autonomous d2F/dtdx.

xdim = data.dim;
ydim = data.ydim;

if ~isempty(p)
    if isempty(data.dfdtdxhan)
        f  = @(x,p) data.dfdthan(p(ydim+1,:), x, p(1:ydim,:), p(ydim+2:end,:));
        ytp = [ y; t ; p ];
        Jtx = coco_ezDFDX('f(x,p)', f, x, ytp);
    else
        [m, n] = size(x);
        Jtx = zeros(m,m,n);
        dfdtdxhan = data.dfdtdxhan;
        for i=1:n
            Jtx(:,:,i) = dfdtdxhan(t(i), x(:,i), y(:,i), p(:,i));
        end
    end
  
else
    if isempty(data.dfdtdxhan)
        f  = @(x,p) data.dfdthan(p(ydim+1,:), x, p(1:ydim,:));
        yt = [ y; t ];
        Jtx = coco_ezDFDX('f(x,p)', f, x, yt);
    else
        [m, n] = size(x);
        Jtx = zeros(m,m,n);
        dfdtdxhan = data.dfdtdxhan;
        for i=1:n
            Jtx(:,:,i) = dfdtdxhan(t(i), x(:,i), y(:,i));
        end
    end

end
end

function Jty = het_DFDTDY(data, t, x, y, p)
%DFDTDX   Vectorized evaluation of non-autonomous d2F/dtdx.

xdim = data.dim;
ydim = data.ydim;

if ~isempty(p)
    if isempty(data.dfdtdyhan)
        f  = @(x,p) data.dfdthan(p(xdim+1,:), p(1:xdim,:), x, p(xdim+2:end,:));
        xtp = [ x; t ; p ];
        Jty = coco_ezDFDX('f(x,p)', f, y, xtp);
    else
        [m, n] = size(y);
        Jty = zeros(size(x,1),m,n);
        dfdtdyhan = data.dfdtdyhan;
        for i=1:n
            Jty(:,:,i) = dfdtdyhan(t(i), x(:,i), y(:,i), p(:,i));
        end
    end
else
    
    if isempty(data.dfdtdyhan)
        f  = @(x,p) data.dfdthan(p(xdim+1,:), p(1:xdim,:), x );
        xt = [ x; t ];
        Jty = coco_ezDFDX('f(x,p)', f, y, xt);
    else
        [m, n] = size(y);
        Jty = zeros(size(x,1),m,n);
        dfdtdyhan = data.dfdtdyhan;
        for i=1:n
            Jty(:,:,i) = dfdtdyhan(t(i), x(:,i), y(:,i) );
        end
    end
  
end

end

function Jtp = het_DFDTDP(data, t, x, y, p)
%DFDTDP   Vectorized evaluation of non-autonomous d2F/dtdp.

xdim = data.dim;
ydim = data.ydim;

  if isempty(data.dfdtdphan)
    f  = @(x,p) data.dfdthan(x(end,:), x(1:xdim,:), x(xdim+1:xdim+ydim,:), p);
    xyt = [ x; y; t ];
    Jtp = coco_ezDFDP('f(x,p)', f, xyt, p);
  else
    [m, n] = size(x);
    o  = size(p,1);
    Jtp = zeros(m,o,n);
    dfdtdphan = data.dfdtdphan;
    for i=1:n
      Jtp(:,:,i) = dfdtdphan(t(i), x(:,i), y(:,i), p(:,i));
    end
  end

end

function Jtt = het_DFDTDT(data, t, x, y, p)
%DFDTDT   Vectorized evaluation of non-autonomous d2F/dt2.
xdim = data.dim;
ydim = data.ydim;

if ~isempty(p)
    if isempty(data.dfdtdthan)
        [m, n] = size(x);
        f  = @(x,p) data.dfdthan(x(1,:), p(1:m,:), p(m+1:m+ydim,:), p(m+ydim+1:end,:));
        xyp = [ x ; y; p ];
        Jtt = reshape(coco_ezDFDX('f(x,p)', f, t, xyp), [m n]);
    else
        [m, n] = size(x);
        Jtt = zeros(m,n);
        dfdtdthan = data.dfdtdthan;
        for i=1:n
            Jtt(:,i) = dfdtdthan(t(i), x(:,i), y(:,i), p(:,i));
        end
    end
else
    
    if isempty(data.dfdtdthan)
        [m, n] = size(x);
        f  = @(x,p) data.dfdthan(x(1,:), p(1:m,:), p(m+1:m+ydim,:));
        xy = [ x ; y];
        Jtt = reshape(coco_ezDFDX('f(x,p)', f, t, xy), [m n]);
    else
        [m, n] = size(x);
        Jtt = zeros(m,n);
        dfdtdthan = data.dfdtdthan;
        for i=1:n
            Jtt(:,i) = dfdtdthan(t(i), x(:,i), y(:,i));
        end
    end
  
end

end





