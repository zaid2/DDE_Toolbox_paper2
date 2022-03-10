function J = duffing_dfdp(t,x,ycn,p)

x1 = x(1,:);
x2 = x(2,:);

y1 = ycn(1,:); 
y2 = ycn(2,:);

mu     = 1;
zeta   = 5*1e-3;

Omega  = p(1,:);
a      = p(2,:);
b      = p(3,:);

J = zeros(2,size(p,1),numel(x1));
J(2,1,:) = -a.*sin(Omega.*t).*t+b.*cos(Omega.*t).*t;
J(2,2,:) =  cos(Omega.*t);
J(2,3,:) =  sin(Omega.*t);



end






