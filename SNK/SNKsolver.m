function [ u,stats ] = SNKsolver(f,Jac,init,PUApprox,varargin)
% SNKsolver
% The Schwarz Newton Krylov  method (SNK) solves nonlinear PDEs by nonlinear
% preconditioning Newton's method via a alternating Schwarz process.
% With this method, the PDE is solved for overlaping subdomains with
% independent grids.
%
%
% INPUT:
%             f: nonlinear residual function f(x,p) for solution x and local
%             approximation p. f(x,p) evaluates the residual on the domain
%             of p.
%
%            Jac: Jacobian function Jac(x,p) for solution x and local
%             approximation p. Jac(x,p) evaluates the residual on the domain
%             of p.
%
%           init: initial guess for all subdomains.
%
%       PUApprox: PUApprox approximation
%
%          tol_n: [rel_tol, abs_tol] relative and absolute tolerance used
%                 for Newtons method, both with SNK and the local nonlinear
%                 solves.
%
% OUTPUT:
%              u: solution of PDE.
%
%        normres: norm of SNK forward evaluation of each iteration.
%
%       normstep: norm of Newton step for each iteration.
%
%          numgm: number of GMRES iterations be iteration.
%
%       normresf: norm of PDE residual of each iteration.
%
% NOTE u,init is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then x = [u1;u2;v1;v2].

opt = solveroptions(varargin{:});
normres = []; normstep = [];  numgm = []; normresf = [];

u = init;

opt_c = opt;  opt_c.inparallel = false;

restime = 0;  % total time for residual computations
jactime = 0;  % total time for Jacobian applications
resbytes = 0;
jacbytes = 0;

numpoints = length(PUApprox);

% solve for the new value using plain Newton
tol_g = 1e-4;
for k = 1:20
	
	%normresf(k) = norm(NKSresidual(u,PUApprox,f,opt));
	
	if opt.twolevel && k > 0
		[c_sol,Jac_hat] = CoarseCorrect(PUApprox,u,f,Jac,opt_c);
	else
		c_sol = 0;
	end
	
	% evaluate the local corrections/solve local nonlinear problems
	ts = tic;
	if opt.inparallel, pp = gcp;  bs = pp.ticBytes;  end
	[z,luJ] = SNKresidual(u+c_sol,PUApprox,f,Jac,opt);
	restime = restime + toc(ts);
	if opt.inparallel, resbytes = resbytes + pp.tocBytes(bs);  end
	z = z + c_sol;
	
	% normalize to be more like integral 2-norm
	%z = z/sqrt(numpoints);
	
	normres(k) = norm(z);
	fprintf('  Residual norm = %.3e\n',normres(k))
	%normresf(k)
	
	if k==1
		stop_tol = norm(z)*opt.reltol + opt.abstol;
	end
	
	if normres(k) < stop_tol, break, end
	
	%tol_g = max(1e-11,1e-12/normres(k))*norm(u);
	if k > 1
		tol_g = max(min(tol_g,1e-4*(normres(k)/normres(k-1))^2),1e-10);
	end
	
	%solve the newton step
	if opt.twolevel && k > 0
		[FJv,FJv_hat] = ComputeJac(PUApprox,Jac,u,opt_c);
		[s,~,~,~,gmhist] = gmres(@Jacobian,-z,[],tol_g,180);
	else
		[s,~,~,~,gmhist] = gmres(@(x)SNKJacobian(PUApprox,luJ,x,opt),-z,[],tol_g,200);
	end
	normstep(k) = norm(s);  numgm(k) = length(gmhist) - 1;
	
	fprintf('  Newton step took %i GMRES iterations\n\n',numgm(k))
	
	% update
	u = u+s;
end

stats.normres = normres;
stats.normstep = normstep;
stats.numgmres = numgm;
stats.normresf = normresf;
stats.restime = restime;
stats.jactime = jactime;
stats.resbytes = sum(resbytes,1);
stats.jacbytes = sum(jacbytes,1);

	function y = Jacobian(w)
		ts = tic;
		c_w = LinearCoarseCorrect( PUApprox, w,Jac_hat,FJv,FJv_hat,opt_c);
		if opt.inparallel, pp = gcp;  bs = pp.ticBytes;  end
		y = c_w + SNKJacobian(PUApprox,luJ,w+c_w,opt);
		jactime = jactime + toc(ts);
		if opt.inparallel, jacbytes = jacbytes + pp.tocBytes(bs);  end
	end

end


