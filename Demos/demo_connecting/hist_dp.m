function J = hist_dp(t,p)

p1 = p(1,:);
p2 = p(2,:);
alpha = p(3,:);
lambda = p(4,:);
v1     = p(5,:);
v2     = p(6,:);
eps = p(7,:);

v = [v1;v2];

% eps = 1e-3;
eps = 7.107e-4;

J = zeros(2,size(p,1),length(t));

J(:,4,:) = eps*(v.*t).*exp(lambda.*t);
J(1,5,:) = eps*exp(lambda.*t);
J(2,6,:) = eps*exp(lambda.*t);
J(:,7,:) = v.*exp(lambda.*t);

end







