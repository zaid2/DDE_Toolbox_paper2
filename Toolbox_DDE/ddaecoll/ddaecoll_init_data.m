function data = ddaecoll_init_data(data, x0, y0, p0)
%COLL_INIT_DATA   Initialize toolbox data for an instance of 'coll'.
%
% Populate remaining fields of the toolbox data structure used by 'coll'
% function objects.
%
% DATA = COLL_INIT_DATA(DATA, X0, P0)
%
% DATA - Toolbox data structure.
% X0   - Initial solution guess for discretized trajectory.
% P0   - Initial solution guess for problem parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_init_data.m 2839 2015-03-05 17:09:01Z fschild $

NTST = data.ddaecoll.NTST; % Number of mesh intervals
NCOL = data.ddaecoll.NCOL; % Degree of polynomial interpolants
dim  = size(x0, 2);    % State-space dimension
ydim = size(y0,2);

pdim = numel(p0);      % Number of problem parameters

data.dim  = dim;       % State-space dimension
data.ydim = ydim;      % Dimension of y variables
data.pdim = pdim;      % Number of problem parameters

bpnum  = NCOL+1;            % Number of basepoints per interval
bpdim  = dim*(NCOL+1);      % Number of basepoint values per interval
bpydim = ydim*(NCOL+1);
xbpnum = (NCOL+1)*NTST;     % Number of basepoints
xbpdim = dim*(NCOL+1)*NTST; % Number of basepoint values
ybpdim = ydim*(NCOL+1)*NTST; % Number of basepoint values
cndim  = dim*NCOL;          % Number of collocation node values per interval

xcnnum = NCOL*NTST;         % Number of collocation nodes
xcndim = dim*NCOL*NTST;     % Number of collocation conditions
ycndim = ydim*NCOL*NTST;    % Number of collocation conditions for y 
cntnum = NTST-1;            % Number of internal boundaries
cntdim = dim*(NTST-1);      % Number of continuity conditions
cnydim = ydim*NCOL;

data.xbpnum  = xbpnum;
data.xcnnum  = xcnnum;
data.xbpdim  = xbpdim;
data.xcndim  = xcndim;
data.ybpdim  = ybpdim;
data.ycndim  = ycndim;
data.xbp_idx = (1:xbpdim)'; % Index array for basepoint values
data.ybp_idx = xbpdim+(1:ybpdim)';
data.T0_idx  = xbpdim+ybpdim+1;
data.T_idx   = xbpdim+ybpdim+2;    % Index for interval length
data.p_idx   = xbpdim+ybpdim+2+(1:pdim)'; % Index array for problem parameters

data.x_shp   = [dim xcnnum]; % Shape for vectorization
data.xbp_shp = [dim xbpnum]; % Shape for vectorization
data.y_shp   = [ydim xcnnum];% Shape for vectorization
data.ybp_shp = [ydim xbpnum]; % Shape for vectorization
data.p_rep   = [1 xcnnum];   % Shape for vectorization
data.pbp_rep   = [1 xbpnum];   % Shape for vectorization

tm = linspace(-1, 1, bpnum)';  
t  = repmat((0.5/NTST)*(tm+1), [1 NTST]);
t  = t + repmat((0:cntnum)/NTST, [bpnum 1]);
data.tbp = t(:)/t(end);

[tc wts] = coll_nodes(NCOL);
wts      = repmat(wts, [dim NTST]);
data.wts    = wts;
data.wts1   = wts(1,:);                          
data.wts2   = spdiags(wts(:), 0, xcndim, xcndim);


data.tm = tm; 
data.tc = tc;

taucn      = kron((0:NTST-1)',ones(NCOL,1))+repmat(0.5*(tc+1),[NTST,1]);
data.taucn     = taucn/NTST;

data.x0_idx = (1:dim)';            % Index array for trajectory end point at t=0
data.x1_idx = xbpdim-dim+(1:dim)'; % Index array for trajectory end point at t=1
data.tbp_idx = setdiff(1:xbpnum, 1+bpnum*(1:cntnum))'; % Index array without duplication of internal boundaries
data.taL_idx = setdiff(1:xcnnum, 1+NCOL*(1:cntnum))';

L   = coll_L(tm, tc); Lp  = coll_Lp(tm, tc);

rows        = reshape(1:xcndim, [cndim NTST]);
rows        = repmat(rows, [bpdim 1]);
cols        = repmat(1:xbpdim, [cndim 1]);
W           = repmat(kron(L, eye(dim)), [1 NTST]);  % Interpolation matrix
Wp          = repmat(kron(Lp, eye(dim)), [1 NTST]); % Interpolation matrix
data.W      = sparse(rows, cols, W);     %===for the segment
data.Wp     = sparse(rows, cols, Wp);

rowsy         = reshape(1:ycndim, [cnydim NTST]);
rowsy        = repmat(rowsy, [bpydim 1]);
colsy        = repmat(1:ybpdim, [cnydim 1]);
Wy            = repmat(kron(L, eye(ydim)), [1 NTST]);  % Interpolation matrix
Wyp           = repmat(kron(Lp, eye(ydim)), [1 NTST]); % Interpolation matrix
data.Wy       = sparse(rowsy,colsy,Wy);
data.Wyp      = sparse(rowsy,colsy,Wyp);

data.dxrows = repmat(reshape(1:xcndim, [dim xcnnum]), [dim 1]);  % Index array for vectorization
data.dxcols = repmat(1:xcndim, [dim 1]);                         % Index array for vectorization
data.dprows = repmat(reshape(1:xcndim, [dim xcnnum]), [pdim 1]); % Index array for vectorization
data.dpcols = repmat(1:pdim, [dim xcnnum]);                      % Index array for vectorization
data.dyrows = repmat(reshape(1:xcndim, [dim xcnnum]), [ydim 1]);  % Index array for vectorization
data.dycols = repmat(1:ycndim, [dim 1]);                          % Index array for vectorization
data.dywdbprows = repmat(reshape(1:(NTST*NCOL+1)*dim, [dim NTST*NCOL+1]), [ydim 1]);  % Index array for vectorization
data.dywdbpcols = repmat(1:(NTST*NCOL+1)*ydim, [dim 1]);                          % Index array for vectorization
data.dtrows = reshape(1:xcndim, [dim xcnnum]);                    % Index array for vectorization
data.dtcols = repmat(1:xcnnum, [dim 1]);                          % Index array for vectorization

temp        = reshape(1:xbpdim, [bpdim NTST]);
Qrows       = [1:cntdim 1:cntdim];
Qcols       = [temp(1:dim, 2:end) temp(cndim+1:end, 1:end-1)];
Qvals       = [ones(cntdim,1) -ones(cntdim,1)];
data.Q      = sparse(Qrows, Qcols, Qvals, cntdim, xbpdim); % Jacobian of continuity conditions
data.dyTpcnt= sparse(cntdim, ybpdim+2+pdim);

end

function [nds wts] = coll_nodes(m)
%COLL_NODES   Compute collocation nodes and integration weights.
%
% Uses eigenvalues and eigenvectors of Jacobi matrix.
%
% [NDS WTS] = COLL_NODES(M)
%
% NDS - Collocation nodes.
% WTS - Quadrature weights.
% M   - Polynomial degree.

n = (1:m-1)';
g = n.*sqrt(1./(4*n.^2-1));
J = -diag(g,1)-diag(g,-1);

[w x] = eig(J);
nds   = diag(x);
wts   = 2*w(1,:).^2;

end

function A = coll_L(ts, tz)
%COLL_L   Evaluation of Lagrange polynomials.
%
% Use high-dimensional arrays for vectorized evaluation.
%
% A = COLL_L(TS, TZ)
% 
% A  - Array of interpolated values.
% TS - Array of basepoints.
% TZ - Array of interpolation points.

q = numel(ts);
p = numel(tz);

zi = repmat(reshape(tz, [p 1 1]), [1 q q]);
sj = repmat(reshape(ts, [1 q 1]), [p 1 q]);
sk = repmat(reshape(ts, [1 1 q]), [p q 1]);

t1 = zi-sk;
t2 = sj-sk;
idx = find(abs(t2)<=eps);
t1(idx) = 1;
t2(idx) = 1;

A = prod(t1./t2, 3);

end

function A = coll_Lp(ts, tz)
%COLL_LP   Evaluation of derivative of Lagrange polynomials.
%
% Use high-dimensional arrays for vectorized evaluation.
%
% A = COLL_LP(TS, TZ)
% 
% A  - Array of interpolated values
% TS - Array of basepoints
% TZ - Array of interpolation points

q = numel(ts);
p = numel(tz);

zi = repmat(reshape(tz, [p 1 1 1]), [1 q q q]);
sj = repmat(reshape(ts, [1 q 1 1]), [p 1 q q]);
sk = repmat(reshape(ts, [1 1 q 1]), [p q 1 q]);
sl = repmat(reshape(ts, [1 1 1 q]), [p q q 1]);

t3 = sj(:,:,:,1)-sk(:,:,:,1);
t4 = zi-sl;
t5 = sj-sl;

idx1 = find(abs(t5)<=eps);
idx2 = find(abs(t3)<=eps);
idx3 = find(abs(sk-sl)<=eps);
t5(union(idx1, idx3)) = 1;
t4(union(idx1, idx3)) = 1;
t3(idx2) = 1;
t3       = 1.0./t3;
t3(idx2) = 0;

A = sum(t3.*prod(t4./t5, 4), 3);

end



