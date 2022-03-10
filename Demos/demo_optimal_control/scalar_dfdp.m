function J = scalar_dfdp(t,x,y,p)

x  = x(1,:);
xd = y(1,:);
u  = y(2,:);
% p  = p(1,:);

J = zeros(1,size(p,1),length(t));
% J(1,1,:) = y(1,:);

% ydot = t.*x+p.*xd+u;
end