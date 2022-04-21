function [data,y]=sys_obj(prob,data,u)

% Copyright (C) Zaid Ahsan

NTST = data.ddaecoll.NTST;
W = data.W; wts = data.wts;
xbp = u(data.xbp_idx);
ubp = u(data.xbp_idx(end)+1:end,1);

xcn = W*xbp; ucn = W*ubp;
gcn = xcn.^2 + ucn.^2;
y = (1/NTST)*wts*gcn(:);


end
