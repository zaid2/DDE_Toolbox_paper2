function J = ddetorus_dfdydy(t, x, y, p)
%TORUS_DFDX   'coll'-compatible encoding of Jacobian of vector field w.r.t. problem variables
% Copyright (C) Zaid Ahsan
x1 = x(1,:);
x2 = x(2,:);

y1 = y(1,:);
y2 = y(2,:);

om1  = p(1,:);
Texc = p(2,:);

om2  = 2*pi./Texc;

r = sqrt(x1.^2+x2.^2);

J = zeros(2,2,2,numel(x1));



end
