function J = ddetorus_dfdtdp(t, x, y, p)
%TORUS_DFDT   'coll'-compatible encoding of vector field.

x1 = x(1,:);
x2 = x(2,:);

y1 = y(1,:);
y2 = y(2,:);

Texc = p(2,:);
om2 =  2*pi./Texc;

r = sqrt(x1.^2+x2.^2);
J = zeros(2,3,numel(x1));
J(1,2,:) = -y1.*r.*(-2*pi./Texc.^2).*sin(om2.*t)-y1.*r.*om2.*cos(om2.*t).*(-2*pi./Texc.^2).*t;
J(2,2,:) = -y2.*r.*(-2*pi./Texc.^2).*sin(om2.*t)-y2.*r.*om2.*cos(om2.*t).*(-2*pi./Texc.^2).*t;

end
