function J = mackey_dfdtdy(t,x,y,p)

x  = x(1,:);
y  = y(1,:);
a  = p(1,:);
b  = p(2,:);

J = zeros(1,1,length(t));
% ydot = a.*y./(1+y.^b)-x;
end




