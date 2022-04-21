function J = scalar_dfdy(t,x,y,p)

% Copyright (C) Zaid Ahsan

x  = x(1,:);
xd = y(1,:);
u  = y(2,:);
p  = p(1,:);

J = zeros(1,2, length(t));

J(1,1,:) = 1; J(1,2,:) = 1;

% ydot = t.*x+p.*xd+u;

end
