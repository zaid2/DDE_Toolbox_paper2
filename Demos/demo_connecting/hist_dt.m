function J = hist_dt(t,p)
% Copyright (C) Zaid Ahsan
p1 = p(1,:);
p2 = p(2,:);
alpha = p(3,:);
lambda = p(4,:);
v1     = p(5,:);
v2     = p(6,:);
eps    = p(7,:);

v = [v1;v2];

% eps = 1e-3;
% eps = 7.107e-4;

J = (eps.*v.*lambda).*exp(lambda.*t);

end
