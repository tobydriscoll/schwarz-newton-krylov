function [ J ] = BurgersJacobian(u,Approx,nu)
%   PDE Jacobian for steady Burgers equation

degs = Approx.degs;
u = reshape(u,degs);

% Locate boundary points
[~,~,~,~,border,~] = FindBoundaryIndex2DSides(degs,Approx.domain,Approx.outerbox);

dx = diffmat(degs(1),1,Approx.domain(1,:));
dy = diffmat(degs(2),1,Approx.domain(2,:));

Ix = eye(degs(1));
Iy = eye(degs(2));
I = eye(prod(degs));

Dx = kron(Iy,dx);
Dy = kron(dy,Ix);
Dxx = kron(Iy,dx^2);
Dyy = kron(dy^2,Ix);

ux = dx*u; uy = u*dy';
u = u(:); ux = ux(:); uy = uy(:);
J = nu*(Dxx+Dyy) - u(:).*(Dx+Dy);
J = J - diag(ux+uy);

% Dirichlet condition
J(border,:) = I(border,:);

end