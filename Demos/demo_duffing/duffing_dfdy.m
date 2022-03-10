function J = duffing_dfdy(t,x,ycn,p)

x1 = x(1,:);
x2 = x(2,:);

y1 = ycn(1,:); 
y2 = ycn(2,:);

mu     = 1;
zeta   = 5*1e-3;

Omega  = p(1,:);
a      = p(2,:);
b      = p(3,:);
gamma  = -0.01;
% gamma = 0;

J = zeros(2,2,numel(x1));
J(2,1,:)= gamma;


end