% CoarseCorrect
% This method solves for the nonlinear coarse correction used in the 2
% level SNK method.
% 
% INPUT:
%      PUApprox: PUApprox approximation
%             v: given solution
%
%             f: nonlinear residual function f(x,p) for solution x and local
%             approximation p. f(x,p) evaluates the residual on the domain
%             of p.
%
%            Jac: Jacobian function Jac(x,p) for solution x and local
%             approximation p. Jac(x,p) evaluates the residual on the domain
%             of p.
%
% OUTPUT:
%          r_er: correction of solution
%          J_v_pls_er: jacobian of J_hat(v+er)
%
% NOTE sol is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then sol = [u1;u2;v1;v2].
function [ r_er,J_v_pls_er ] = CoarseCorrect( PUApprox, v,f,Jac,opt)

gd = opt.griddepth;
tol_c = opt.coarsetol;
numpoints = length(PUApprox);

v_hat = [];
for i=1:opt.numcomp
    v_hat = [v_hat;PUApprox.Fine2Coarse(v((1:numpoints) + (i-1)*numpoints),gd)];
end

r = NKSresidual(v,PUApprox,f,opt);

r_hat = [];
for i=1:opt.numcomp
    r_hat = [r_hat;PUApprox.Fine2Coarse(r((1:numpoints) + (i-1)*numpoints),gd)];
end

PUApprox.Coarsen();

r_hat = r_hat - NKSresidual(v_hat,PUApprox,f,opt);

RES = @(er)Residual(er,v_hat,r_hat,PUApprox,f,opt);
JAC = @(er)CoarseASJac(PUApprox,Jac,er,v_hat);

fprintf('  Coarse solution...')
options = optimoptions(@fsolve,...
    'algorithm','trust-region-dogleg','derivativecheck','off','SpecifyObjectiveGradient',true,'MaxIterations',200,'FunctionTolerance',tol_c,'Display','final');
[er,~,flag] = fsolve(@(er)sol_and_jac(@(er)RES(er),@(er)JAC(er),er),zeros(size(v_hat)),options);
er = er(:,end);
if flag <=0
    warning('fsolve exited with flag %i',flag)
end
fprintf('finished.\n')

J_v_pls_er = JAC(er);

r_er = [];

numcoarse = numel(r_hat)/opt.numcomp;
for i=1:opt.numcomp
    r_er = [r_er;PUApprox.Coarse2Fine(er((1:numcoarse) + (i-1)*numcoarse))];
end

PUApprox.Refine();

    function [ sol,J ] = sol_and_jac( f,jac,u )
        sol = f(u);
        J = jac(u);
    end

end

function F = Residual(er,v,r,PUApprox,evalf,opt)
    F = NKSresidual(v+er,PUApprox,evalf,opt);
    F = F + r;
end

