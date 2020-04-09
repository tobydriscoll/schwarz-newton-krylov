function F = Burgers(u,Approx,nu,bound_f)
%   PDE residual for steady Burgers equation

degs = Approx.degs;
u = reshape(u,degs);

% Locate boundary points
[~,~,~,~,border,~] = FindBoundaryIndex2DSides(degs,Approx.domain,Approx.outerbox);

Dx = diffmat(degs(1),1,Approx.domain(1,:));
Dy = diffmat(degs(2),1,Approx.domain(2,:));

% PDE residual computed here
ux = Dx*u; uxx = Dx*ux; uy = u*Dy'; uyy = uy*Dy';
F = nu*(uxx+uyy) - u.*(ux+uy);

P = Approx.points();

% Dirichlet condition applied here
F(border) = u(border) - bound_f(P(border,1),P(border,2)); 

F = F(:);

end
