function y = hist(t,p)

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

y = (eps.*v).*exp(lambda.*t);


end


% eta = p(1,1);
% A0=[0 1;1 eta];
% [U,V]=eig(A0);
% 
% vkminus = U(:,2);
% lambdaminus = V(2,2);

% vkminus = [-0.858642953046343; -0.512574169446585];
% lambdaminus = 0.596958453601750;
% hold on 
% plot(t,y)
