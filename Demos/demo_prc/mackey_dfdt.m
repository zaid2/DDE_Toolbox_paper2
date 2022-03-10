function J = mackey_dfdt(t,x,y,p)

x  = x(1,:);
y  = y(1,:);
a  = p(1,:);
b  = p(2,:);

J = zeros(1,length(t));
% ydot = a.*xd./(1+xd.^b)-x;
end