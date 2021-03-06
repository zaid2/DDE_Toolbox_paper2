function J = mackey_dfdydy(t,x,y,p)
% Copyright (C) Zaid Ahsan
x  = x(1,:);
y  = y(1,:);
a  = p(1,:);
b  = p(2,:);

J = zeros(1,1,1, length(t));

J(1,1,1,:) = a.*(-1./(1+y.^b).^2).*(b.*(y.^(b-1)))+...
             a.*(-1./(1+y.^b).^2).*(b.*(y.^(b-1)))+...
             a.*y.*(2./(1+y.^b).^3).*(b.*(y.^(b-1))).*(b.*(y.^(b-1)))+...
             a.*y.*(-1./(1+y.^b).^2).*(b.*(b-1).*y.^(b-2));
         
% ydot = a.*y./(1+y.^b)-x;
end
