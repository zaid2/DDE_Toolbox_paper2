function J = mackey_dfdx(t,x,y,p)

x  = x(1,:);
y  = y(1,:);
a  = p(1,:);
b  = p(2,:);

J = -1*ones(1,1, length(t));

% ydot = a.*y./(1+y.^b)-x;

end