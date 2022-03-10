function [data, J] = glue_coup_du(prob, data, u)

alpha     = u(end-1);
T         = u(end);

J = zeros(data.nsegs,data.nsegs+2);

for i=1:data.nsegs

J(i,i) = 1;
J(i,data.nsegs+1) = -1/T;
J(i,data.nsegs+2) = alpha/T^2;

end




% dJ_incorr = full(J);
% [data, dJ_corr] = coco_ezDFDX('f(o,d,x)', prob, data, @glue_coup, u);

end









