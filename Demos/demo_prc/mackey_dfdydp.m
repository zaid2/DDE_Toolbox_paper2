function J = mackey_dfdydp(t,x,y,p)
% Copyright (C) Zaid Ahsan
x  = x(1,:);
y  = y(1,:);
a  = p(1,:);
b  = p(2,:);

J = zeros(1,1,3, length(t));
J(1,1,1,:) = 1./(1+y.^b)+y.*(-1./(1+y.^b).^2).*(b.*(y.^(b-1)));

J(1,1,2,:) = a.*(-1./(1+y.^b).^2).*(log(y).*(y.^b))+...
             a.*y.*((2./(1+y.^b).^3).*(y.^b).*log(y)).*(b.*(y.^(b-1)))+...
             a.*y.*(-1./(1+y.^b).^2).*(y.^(b-1))+...
             a.*y.*(-1./(1+y.^b).^2).*b.*(y.^(b-1)).*log(y);
         
         
% ydot = a.*y./(1+y.^b)-x;

end

