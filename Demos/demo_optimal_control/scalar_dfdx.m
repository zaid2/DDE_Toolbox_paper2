function J = scalar_dfdx(t,x,y,p)

x  = x(1,:);
xd = y(1,:);
u  = y(2,:);
% p  = p(1,:);

J = zeros(1,1, length(t));

J(1,1,:) = t;

% ydot = t.*x+p.*xd+u;

end