addpath ../PUChebfun
addpath ../SNK

domain = [-1 1;-1 1];
cheb_struct.domain = domain;
cheb_struct.degs = [33 33];
cheb_struct.cdegs = [9 9];
cheb_struct.split_flag = [true true];
cheb_struct.tol = 1e-4;

% Solution tree setup
Tree = ChebPatch(cheb_struct);
Tree = Tree.split(1);
for dim = [2 1 2]
	Tree.split(dim);
end
numleaves = length(collectLeaves(Tree))

F = PUchebfun(Tree);
F.sample(@(x,y) zeros(size(x)));

setInterpMatrices(F,true);
bound_f = @(x,y) atan((cos(pi*3/16)*x+sin(pi*3/16)*y)*1);
nu = 1/100;
f = @(u,leaf) Burgers(u,leaf,nu,bound_f);
Jac = @(u,leaf) BurgersJacobian(u,leaf,nu);
F.Setvalues(bound_f);
init = F.Getvalues();

opt = [];
opt.reltol = 1e-10;
opt.abstol = 1e-11;
opt.inparallel = false;
opt.numwork = 2;   % number of workers
opt.coarsetol = 1e-5;
opt.griddepth = 2;
opt.twolevel = true;  % only for SNK
opt.solver = @NKSsolver;
%opt.solver = @SNKsolver;
opt = solveroptions(opt)
