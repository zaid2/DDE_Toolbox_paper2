function J = mackey_dfdpdp(t,x,y,p)

% Copyright (C) Zaid Ahsan

x  = x(1,:);
y  = y(1,:);
a  = p(1,:);
b  = p(2,:);

J = zeros(1,3,3,length(t));

J(1,2,1,:) = y.*(-1./(1+y.^b).^2).*log(y).*(y.^b);
J(1,1,2,:) = y.*(-1./(1+y.^b).^2).*((y.^b).*log(y));
J(1,2,2,:) = a.*y.*(2./(1+y.^b).^3).*(log(y).*(y.^b)).*log(y).*(y.^b)+...
             a.*y.*(-1./(1+y.^b).^2).*log(y).*(log(y).*(y.^b));
end






