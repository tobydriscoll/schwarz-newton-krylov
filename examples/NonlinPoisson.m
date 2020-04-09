function F = NonlinPoisson(u,Approx)
%   PDE residual for a nonlinear Poisson equation

degs = Approx.degs;
u = reshape(u,degs);

west = false(degs); west(1,2:end-1) = true;
east = false(degs); east(end,2:end-1) = true;
south = false(degs); south(:,1) = true;
north = false(degs); north(:,end) = true;

Dx = diffmat(degs(1),1,Approx.domain(1,:));
Dy = diffmat(degs(2),1,Approx.domain(2,:));

ux = Dx*u; uy = u*Dy';
[X,Y] = ndgrid(Approx.leafGrids{:});

F = X.*sin(pi*Y) - Dx*((1+u.^2).*ux) - ((1+u.^2).*uy)*Dy';

F(east) = u(east)-1;  % Dirichet
F(west) = ux(west);   % Neumann
F(south | north) = uy(south | north);  % Neumann

F = F(:);
end
