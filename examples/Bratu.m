function F = Bratu(u,Approx,lambda)
%   PDE residual for Bratu equation

degs = Approx.degs;
u = reshape(u,degs);

% Locate boundary points
[~,~,~,~,border,~] = FindBoundaryIndex2DSides(degs,Approx.domain,Approx.outerbox);

Dx = diffmat(degs(1),1,Approx.domain(1,:));
Dy = diffmat(degs(2),1,Approx.domain(2,:));

% PDE residual computed here
uxx = Dx^2*u; uyy = u*(Dy')^2;
F = uxx + uyy + lambda*exp(u);

F(border) = u(border);  % Dirichlet condition
F = F(:);

end

