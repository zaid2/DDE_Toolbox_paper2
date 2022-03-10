function J = scalar_dfdt(t,x,y,p)

x  = x(1,:);
xd = y(1,:);
u  = y(2,:);
% p  = p(1,:);



J = zeros(1,length(t));
J(1,:) = x;

% ydot = t.*x+p.*xd+u;
end