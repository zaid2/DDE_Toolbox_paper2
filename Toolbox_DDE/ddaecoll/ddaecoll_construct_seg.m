function prob = ddaecoll_construct_seg(prob, tbid, data, sol)
%COLL_CONSTRUCT_SEG   Append an instance of 'coll' to problem.
%
% Add collocation and continuity conditions, monitor functions that
% evaluate to the problem parameters, and corresponding inactive
% continuation parameters.
%
% PROB = COLL_CONSTRUCT_SEG(PROB, TBID, DATA, SOL)
%
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.
% SOL  - Initial solution guess.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_construct_seg.m 2839 2015-03-05 17:09:01Z fschild $

prob = coco_add_func(prob, tbid, @ddaecoll_F, @ddaecoll_DFDU, data, 'zero', ...
  'u0', sol.u);

% prob = coco_add_func(prob, tbid, @ddaecoll_F,  data, 'zero', ...
%   'u0', sol.u);

if ~isempty(data.pnames) % Optional monitor functions
  fid  = coco_get_id(tbid, 'pars');
  uidx = coco_get_func_data(prob, tbid, 'uidx'); % Context-dependent index array
  prob = coco_add_pars(prob, fid, uidx(data.p_idx), data.pnames);
end
prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');

end

function [data y] = ddaecoll_F(prob, data, u)
%COLL_F   Collocation zero problem.
%
% Expects vectorized encoding of vector field as function of two arguments.

x  = u(data.xbp_idx); % Extract basepoint values
y  = u(data.ybp_idx);
T0 = u(data.T0_idx);
T  = u(data.T_idx);   % Extract interval length
p  = u(data.p_idx);   % Extract problem parameters

NTST = data.ddaecoll.NTST;


xx  = reshape(data.W*x, data.x_shp); % Values at collocation nodes
yy  = reshape(data.Wy*y, data.y_shp); 
tcn = T*data.taucn+T0;

if data.pdim~=0
pp  = repmat(p, data.p_rep);
ode = data.fhan(tcn',xx, yy, pp);
else
ode = data.fhan(tcn',xx, yy);    
end
ode = 2*NTST*data.Wp*x-T*ode(:);               % Collocation conditions
cnt = data.Q*x;                                % Continuity conditions

y = [ode; cnt];

end

function [data J] = ddaecoll_DFDU(prob, data, u)
%COLL_DFDU   Linearization of collocation zero problem.
%
% Expects vectorized encoding of Jacobians of vector field with respect to
% its arguments. When dfdxhan and/or dfdphan are empty, approximate
% Jacobians are obtained using numerical differentiation.

x  = u(data.xbp_idx); % Extract basepoint values
y  = u(data.ybp_idx);
T0 = u(data.T0_idx);
T  = u(data.T_idx);   % Extract interval length

xx = reshape(data.W*x, data.x_shp); % Values at collocation nodes
yy = reshape(data.Wy*y, data.y_shp); 

if data.pdim~=0
p  = u(data.p_idx);   % Extract problem parameters
pp = repmat(p, data.p_rep);
else
pp = [];    
end

tcn  = data.taucn*T+T0;

dxode = het_DFDX(data,tcn,xx,yy,pp);
dxode = sparse(data.dxrows, data.dxcols, dxode(:));
dxode = 2*data.ddaecoll.NTST*data.Wp-T*dxode*data.W; % W.r.t. basepoint values

if ~isempty(yy)
dyode = het_DFDY(data,tcn,xx,yy,pp);
dyode = sparse(data.dyrows, data.dycols, dyode(:));
dyode = -T*dyode*data.Wy;
else
    dyode = [];
end

dtode  = het_DFDT(data,tcn,xx,yy,pp);
dtode  = sparse(data.dtrows, data.dtcols, dtode(:));
dt0ode = T*dtode;
dT0ode = -dt0ode*ones(data.xcnnum,1);

if ~isempty(pp)
dTode = data.fhan(tcn', xx, yy, pp);
else
dTode = data.fhan(tcn', xx, yy);    
end
dTode  = -dTode(:)-dt0ode*data.taucn; 


if ~isempty(pp)
dpode = het_DFDP(data,tcn,xx,yy,pp);
dpode = sparse(data.dprows, data.dpcols, dpode(:));
dpode = -T*dpode; % W.r.t. problem parameters
else
dpode = [];    
end
J = [dxode dyode dT0ode dTode dpode; data.Q data.dyTpcnt]; % data.dyTpcnt = [0...0]

% J_incorr = full(J);
% [data, J_corr] = coco_ezDFDX('f(o,d,x)', prob, data, @ddaecoll_F, u);
% check = J_corr-J_incorr;
% max(abs(check(:)))

end

function Jx = het_DFDX(data,tcn,xx,yy,pp)
ydim = data.ydim;
dim  = data.dim;

if ~isempty(pp)
    if isempty(data.dfdxhan)
        f  = @(x,p) data.fhan(p(ydim+1,:), x, p(1:ydim,:), p(ydim+2:end,:));
        ytp = [yy; tcn'; pp];
        Jx = coco_ezDFDX('f(x,p)v', f, xx, ytp);
    else
        Jx = data.dfdxhan(tcn', xx, yy, pp);
    end

else
  
    if isempty(data.dfdxhan)
        f  = @(x,p) data.fhan(p(ydim+1,:), x, p(1:ydim,:));
        yt = [yy; tcn'];
        Jx = coco_ezDFDX('f(x,p)v', f, xx, yt);
    else
        Jx = data.dfdxhan(tcn', xx, yy);
    end

end

end


function Jy = het_DFDY(data,tcn,xx,yy,pp)
ydim = data.ydim;
dim  = data.dim;

if ~isempty(pp)
    if isempty(data.dfdyhan)
        f  = @(x,p) data.fhan(p(dim+1,:), p(1:dim,:), x, p(dim+2:end,:));
        xtp = [xx; tcn'; pp];
        Jy = coco_ezDFDX('f(x,p)v', f, yy, xtp);
    else
        Jy = data.dfdyhan(tcn', xx, yy, pp);
    end

else
  
    if isempty(data.dfdyhan)
        f  = @(x,p) data.fhan(p(dim+1,:), p(1:dim,:), x);
        xt = [xx; tcn'];
        Jy = coco_ezDFDX('f(x,p)v', f, yy, xt);
    else
        Jy = data.dfdyhan(tcn', xx, yy);
    end

end

end

function Jt = het_DFDT(data,tcn,xx,yy,pp)
ydim = data.ydim;
dim  = data.dim;

if ~isempty(pp)
    if isempty(data.dfdthan)
        f  = @(x,p) data.fhan(x, p(1:dim,:), p(dim+1:dim+ydim,:), p(dim+ydim+1:end,:));
        xyp = [xx; yy; pp];
        Jt = coco_ezDFDX('f(x,p)v', f, tcn', xyp);
    else
        Jt = data.dfdthan(tcn', xx, yy, pp);
    end

else
  
    if isempty(data.dfdthan)
        f  = @(x,p) data.fhan(x, p(1:dim,:), p(dim+1:dim+ydim,:));
        xy = [xx; yy];
        Jt = coco_ezDFDX('f(x,p)v', f, tcn', xy);
    else
        Jt = data.dfdthan(tcn', xx, yy);
    end

end
end

function Jp = het_DFDP(data,tcn,xx,yy,pp)
ydim = data.ydim;
dim  = data.dim;

if isempty(data.dfdphan)
    f  = @(x,p) data.fhan(p(dim+ydim+1,:), p(1:dim,:), p(dim+1:dim+ydim,:), x);
    xyt = [xx; yy; tcn'];
    Jp = coco_ezDFDX('f(x,p)v', f, pp, xyt);
else
    Jp = data.dfdphan(tcn', xx, yy, pp);
end

end





