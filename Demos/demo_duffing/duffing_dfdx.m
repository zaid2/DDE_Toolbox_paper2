function J = duffing_dfdx(t,x,ycn,p)

x1 = x(1,:);
x2 = x(2,:);

y1 = ycn(1,:); 
y2 = ycn(2,:);

mu     = 1;
zeta   = 5*1e-3;

Omega  = p(1,:);
a      = p(2,:);
b      = p(3,:);

J = zeros(2,2,numel(x1));
J(1,2,:) = 1;
J(2,1,:) = -1-3*mu.*x1.^2;
J(2,2,:) = -2*zeta;




end




