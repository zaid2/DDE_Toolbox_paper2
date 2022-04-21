function y = ddetorus(t, x, y, p)
%TORUS   'coll'-compatible encoding of vector field.
% Copyright (C) Zaid Ahsan
x1 = x(1,:);
x2 = x(2,:);

y1 = y(1,:);
y2 = y(2,:);

om1  = p(1,:);
Texc = p(2,:);

om2 = 2*pi./Texc;

r  = sqrt(x1.^2+x2.^2);
% y(1,:) = -fr.*x2+x1.*(1+r.*(cos(om.*t)-1));
% y(2,:) = fr.*x1+x2.*(1+r.*(cos(om.*t)-1));

y(1,:) = -om1.*x2+y1.*(1+r.*(cos(om2.*t)-1));
y(2,:) = om1.*x1+y2.*(1+r.*(cos(om2.*t)-1));

end
