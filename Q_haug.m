function out = Q_haug(x,y,z)
%Computes the Haugazeau projection onto the intersection of
%half-spaces defined by x, y, and z; details in "Convex
%Analysis and Monotone Operator Theory", 2nd. ed. by Bauschke and
%Combettes, 2017 -- in Definition 29.24.
xi = dot(x(:)-y(:),y(:)-z(:));
mu = norm(x(:)-y(:)).^2;
nu = norm(y(:)-z(:)).^2;
rho = mu*nu - (xi.^2); %always positive
xinu = xi*nu;
if rho==0
    if xi<0
        fprintf('ERROR: C is empty set')
    else
        out = z;
    end
elseif xinu >= rho
    out = x + (1+(xi/nu)).*(z-y);
else %xinu < rho
    out = y + (nu/rho).*( xi*(x-y) + mu.*(z-y) );
end
end
