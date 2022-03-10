function [data,dJ] = sys_obj_du(prob,data,u)

NTST = data.ddaecoll.NTST;
W = data.W; wts = data.wts;
xbp = u(data.xbp_idx);
ubp = u(data.xbp_idx(end)+1:end,1);

xcn = W*xbp; ucn = W*ubp;
% gcn = xcn.^2 + ucn.^2;

gdxcn = 2*diag(xcn(:)); gducn = 2*diag(ucn(:));
dJ = zeros(1,length(u));

dJ(1,1:data.xbp_idx(end)) = (1/NTST)*wts*gdxcn*W;
dJ(1,data.xbp_idx(end)+1:end) = (1/NTST)*wts*gducn*W;

% dJ_incorr = dJ;
% [data, dJ_corr] = coco_ezDFDX('f(o,d,x)', prob, data, @sys_obj, u);
% check = dJ_corr-dJ_incorr;


end