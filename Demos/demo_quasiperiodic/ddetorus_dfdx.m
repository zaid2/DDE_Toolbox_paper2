function J = ddetorus_dfdx(t, x, y, p)
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

  J = zeros(2,2,numel(x1));

  J(1,1,:) = y1.*x1.*(cos(om2.*t)-1)./r;
  J(1,2,:) = -om1 + y1.*(x2.*(cos(om2.*t)-1))./r;
  J(2,1,:) = om1 + y2.*x1.*(cos(om2.*t)-1)./r;
  J(2,2,:) = y2.*x2.*(cos(om2.*t)-1)./r;


end
