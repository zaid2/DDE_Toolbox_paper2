function [data, J] = glue_coup_du(prob, data, u)


gamma1e = u(1);
T       = u(2);
alpha   = u(3);

J = zeros(1,3);

J(1,1)=1; J(1,2) = alpha/T^2; J(1,3)=-1/T;


% [data, dJ_corr] = coco_ezDFDX('f(o,d,x)', prob, data, @glue_coup, u);
% check = dJ_corr-J;


end