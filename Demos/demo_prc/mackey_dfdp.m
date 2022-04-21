function J = mackey_dfdp(t,x,y,p)

% Copyright (C) Zaid Ahsan

x  = x(1,:);
y  = y(1,:);
a  = p(1,:);
b  = p(2,:);

J = zeros(1,3,length(t));
J(1,1,:) = y./(1+y.^b);
J(1,2,:) = a.*y.*(-1./(1+y.^b).^2).*log(y).*(y.^b);

% ydot = a.*xd./(1+xd.^b)-x;
end
