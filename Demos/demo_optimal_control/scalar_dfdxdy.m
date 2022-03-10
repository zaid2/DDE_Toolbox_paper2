function J = scalar_dfdxdy(t,x,y,p)

x  = x(1,:);
xd = y(1,:);
u  = y(2,:);
p  = p(1,:);

J = zeros(1,1,2,numel(x));

% ydot = t.*x+p.*xd+u;
end




