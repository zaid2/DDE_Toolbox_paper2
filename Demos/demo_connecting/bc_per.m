function [data, y] = bc_per(prob, data, u)

% Copyright (C) Zaid Ahsan

T0  = u(1);
x02 = u(2);

Tcr = u(end);

xbp_2dim = u(2:end-1,1);

NTST = data.ddaecoll.NTST;
mesh = 0:1/NTST:1; 
m = data.ddaecoll.NCOL;

L = zeros(1,NTST*(m+1));
k = floor(interp1(mesh,(1:NTST+1),Tcr));
tc = 2*NTST*(Tcr-mesh(k))-1;

L(1,1+(m+1)*(k-1):(m+1)*k) = coll_L(data.tm,tc);

y = [T0; x02; L*xbp_2dim];
% y = [T0; x02];
end
