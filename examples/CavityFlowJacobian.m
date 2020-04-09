function [J] = CavityFlowJacobian(Re,y,leaf)
%   PDE Jacobian for regularized lid-driven cavity flow

degs = leaf.degs;
Len = prod(degs);

% Locate boundary points
[~,out_border,in_border,~,~] = FindBoundaryIndex2DSides(degs,leaf.domain,leaf.outerbox);

border = out_border | in_border;

dx = diffmat(degs(1),1,leaf.domain(1,:));
dy = diffmat(degs(2),1,leaf.domain(2,:));

Ix = eye(degs(1));
Iy = eye(degs(2));

w = zeros(degs);

u(:) = y(1:Len);         % x velocity
v(:) = y(Len+(1:Len));   % y velocity
w(:) = y(2*Len+(1:Len)); % vorticity

wx = dx*w; wy = w*dy';
wx = wx(:); wy = wy(:);

Dx = kron(Iy,dx);  Dxx = kron(Iy,dx^2);
Dy = kron(dy,Ix);  Dyy = kron(dy^2,Ix);

I = eye(Len);
Z = zeros(Len);
J = zeros(3*Len,3*Len);
iu = 1:Len;  iv = Len+iu;  iw = Len+iv;

J(iu,iu) = -(Dxx+Dyy);  J(iu,iw) = -Dy;
J(iu(border),:) = 1*[I(border,:) Z(border,:) Z(border,:)];

J(iv,iv) = -(Dxx+Dyy); J(iv,iw) = Dx;
J(iv(border),:) = 1*[Z(border,:) I(border,:) Z(border,:)];

J(iw,iu) = diag(wx);   J(iw,iv) = diag(wy);  J(iw,iw) = -1/Re*(Dxx+Dyy)+(u.*Dx+v.*Dy);
J(iw(border),:) = 1*[Dy(border,:) -Dx(border,:) I(border,:)];

end

