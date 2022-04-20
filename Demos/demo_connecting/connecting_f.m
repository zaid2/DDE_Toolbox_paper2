function ydot = connecting_f(t,x,y,p)

% Copyright (C) Zaid Ahsan

x1 = x(1,:);
x2 = x(2,:);

y1 = y(1,:); 
y2 = y(2,:);

p1   = p(1,:);
p2   = p(2,:);


ydot = [x2;x1-x1.*y1+p2.*x2+p1.*x1.*x2];


end


