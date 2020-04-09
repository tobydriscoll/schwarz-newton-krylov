function F = CavityFlow(Re,y,leaf,steep)
%   PDE residual for regularized lid-driven cavity flow

degs = leaf.degs;
Len = prod(degs);

% Locate boundary points
[out_border_s,~,~,~,border] = FindBoundaryIndex2DSides(degs,leaf.domain,leaf.outerbox);
north = out_border_s{4};

Dx = diffmat(degs(1),1,leaf.domain(1,:));
Dy = diffmat(degs(2),1,leaf.domain(2,:));

u = zeros(degs);
v = zeros(degs);
w = zeros(degs);

u(:) = y(1:Len);         % x velocity
v(:) = y(Len+(1:Len));   % y velocity
w(:) = y(2*Len+(1:Len)); % vorticity

ux = Dx*u; uxx = Dx*ux; uy = u*Dy'; uyy = uy*Dy';
vx = Dx*v; vxx = Dx*vx; vy = v*Dy'; vyy = vy*Dy';
wx = Dx*w; wxx = Dx*wx; wy = w*Dy'; wyy = wy*Dy';

P = leaf.points();

% First equation
f1 = -(uxx+uyy) - wy;
f1 = f1(:);

% classic boundary regularization
bump = 1 - SideBumpFunc(P(north,1),[0 1],steep) ...
	- SideBumpFunc(1-P(north,1),[0 1],steep);
f1(border) = u(border);
f1(north) = 1*(f1(north) - bump);

% Second equation
f2 = -(vxx+vyy)+wx;
f2 = f2(:);
f2(border) = 1*v(border);

% Third equation
f3 = -1/Re*(wxx+wyy)+(u.*wx+v.*wy);
f3 = f3(:);
f3(border) = 1*(w(border) + uy(border) - vx(border));

F =[f1;f2;f3];


end



