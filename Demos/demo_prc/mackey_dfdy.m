function J = mackey_dfdy(t,x,y,p)

x  = x(1,:);
y  = y(1,:);
a  = p(1,:);
b  = p(2,:);

J = zeros(1,1, length(t));

J(1,1,:) = a.*1./(1+y.^b)+a.*y.*(-1./(1+y.^b).^2).*(b.*(y.^(b-1)));

% ydot = a.*y./(1+y.^b)-x;

end