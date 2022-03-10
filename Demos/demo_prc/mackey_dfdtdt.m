function J = mackey_dfdtdt(t,x,y,p)

x  = x(1,:);
xd = y(1,:);
a  = p(1,:);
b  = p(2,:);

J = zeros(1,length(t));

% ydot = a.*xd./(1+xd.^b)-x;
end




