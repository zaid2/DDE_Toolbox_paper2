function J = connecting_dfdx(t,x,y,p)

% Copyright (C) Zaid Ahsan

x1 = x(1,:);
x2 = x(2,:);

y1 = y(1,:); 
y2 = y(2,:);

p1  = p(1,:);
p2   = p(2,:);

J = zeros(2,2,numel(x1));
J(1,2,:) = 1;
J(2,1,:) = 1-y1+p1.*x2;
J(2,2,:) = p2+p1.*x1;

end




