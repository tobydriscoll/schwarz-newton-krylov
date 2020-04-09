addpath ../PUChebfun
addpath ../SNK

domain = [0 1;0 1];
cheb_struct.domain = domain;
cheb_struct.degs = [33 33];
cheb_struct.cdegs = [9 9];
cheb_struct.split_flag = [true true];
cheb_struct.tol = 1e-4;

F = PUchebfun(@(x,y)exp(-y.^10./(1-y.^10)),[0 1;0 1],'Degree',[33 33],'CoarseDegree',[9 9],'tol',1e-5,'SplitAll',true); F.reset() 
F.sample(@(x,y) zeros(size(x)));

setInterpMatrices(F,true);

Re = 100;
steep = 0.4;  
f = @(u,leaf) CavityFlow(Re,u,leaf,steep);
Jac = @(u,leaf) CavityFlowJacobian(Re,u,leaf);
  
bump = @(x,y) y.*(1 - SideBumpFunc(x,[0 1],steep) - SideBumpFunc(1-x,[0 1],steep));
F.Setvalues(bump); %set values
F.sample(bump); %set coeffs
    
u = F.Getvalues();
Fy = diff(F,2,1);
w = -(Fy.Getvalues());
v = zeros(length(F),1);
init = [u;v;w];

opt = [];
opt.numcomp = 3;  % # of solution components
opt.reltol = 1e-10;
opt.abstol = 1e-11;
opt.inparallel = false;  
opt.numwork = 2;   % number of workers
opt.coarsetol = 1e-4;
opt.griddepth = 2;
opt.twolevel = true;
%opt.solver = @NKSsolver;
opt.solver = @SNKsolver;
opt = solveroptions(opt)