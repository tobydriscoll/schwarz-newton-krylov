function J = NonlinPoissonJac(u,Approx)
%   PDE Jacobian for a nonlinear Poisson equation

degs = Approx.degs;

west = false(degs); west(1,2:end-1) = true;
east = false(degs); east(end,2:end-1) = true;
south = false(degs); south(:,1) = true;
north = false(degs); north(:,end) = true;

dx = diffmat(degs(1),1,Approx.domain(1,:));
dy = diffmat(degs(2),1,Approx.domain(2,:));

Ix = eye(degs(1));  Iy = eye(degs(2));
I = eye(prod(degs));

Dx = kron(Iy,dx);  Dy = kron(dy,Ix);

u = reshape(u,degs);
ux = dx*u; uy = u*dy';

u = u(:); ux = ux(:); uy = uy(:);

J = -Dx*(diag(1+u.^2)*Dx+diag(ux)*(diag(2*u)))... 
	- Dy*(diag(1+u.^2)*Dy+diag(uy)*(diag(2*u)));

J(east,:) = I(east,:);
J(west,:) = Dx(west,:);
J(north | south,:) = Dy(north | south,:);

end