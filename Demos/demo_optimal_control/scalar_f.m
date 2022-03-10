function ydot = scalar_f(t,x,y,p)

x  = x(1,:);
xd = y(1,:);
u  = y(2,:);
% p  = p(1,:);

ydot = t.*x+xd+u;
end