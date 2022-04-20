function y = duffing(t,x,ycn,p)

% Copyright (C) Zaid Ahsan

x1 = x(1,:);
x2 = x(2,:);

y1 = ycn(1,:); 
y2 = ycn(2,:);

mu     = 1;
zeta   = 5*1e-3;

Omega  = p(1,:);
a      = p(2,:);
b      = p(3,:);
gamma = -0.01;
% k = 0;

y = [x2;-2*zeta.*x2-x1-mu.*x1.^3+gamma.*y1+a.*cos(Omega.*t)+b.*sin(Omega.*t)];




end


