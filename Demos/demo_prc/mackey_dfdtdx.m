function J = mackey_dfdtdx(t,x,y,p)
% Copyright (C) Zaid Ahsan
x  = x(1,:);
y  = y(1,:);
a  = p(1,:);
b  = p(2,:);

J = zeros(1,1, length(t));
% ydot = a.*y./(1+y.^b)-x;
end




