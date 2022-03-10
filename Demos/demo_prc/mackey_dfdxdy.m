function J = mackey_dfdxdy(t,x,y,p)

x  = x(1,:);
y  = y(1,:);
a  = p(1,:);
b  = p(2,:);

J = zeros(1,1,1,numel(t));

% ydot = a.*y./(1+y.^b)-x;
end




