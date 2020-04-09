function [z,luJ] = SNKresidual(u,PUApprox,evalF,Jac,opt)
% Residual for the Schwarz Newton Krylov (SNK) method.
%
% INPUT:
%
%            sol: solution used for computing the residual.
%
%            PUApprox: PUApprox approximation
%
%            evalF: nonlinear residual function f(x,p) for solution x and local
%             approximation p. f(x,p) evaluates the residual on the domain
%             of p.
%
%            Jac: Jacobian function f(x,p) for solution x and local
%             approximation p. f(x,p) evaluates the residual on the domain
%             of p.
%
%          tol_n: [rel_tol, abs_tol] relative and absolute tolerance used
%                 for local Newtons method.
% OUTPUT:
%          z: correction of solution
%
%      l,u,p: cell array of LU decomposition and ordering for local
%             Jacobians.
%
%          J: cell array of local Jacobians
%
% NOTE sol is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then sol = [u1;u2;v1;v2].

tol_n = [opt.reltol opt.abstol];

leaf = PUApprox.leafArray;
numcomp = opt.numcomp;
numleaves = length(leaf);
u = reshape(u,length(PUApprox),numcomp);

for k=1:numleaves
	degs{k} = leaf{k}.degs;
	
	%This function returns the logical indicies of the gamma and outer
	%boundry interface. Out put is given for all indicies, as well as the
	%indicies along each of the sides
	[~,~,in_border{k},~] = FindBoundaryIndex2DSides(degs{k},leaf{k}.domain,leaf{k}.outerbox);
	if ~PUApprox.iscoarse
		diff{k} = leaf{k}.Binterp*u;
	else
		diff{k} = leaf{k}.CBinterp*u;
	end
	
end

% Break data into cells for parfor distribution
siz = prod( cat(1,degs{:}), 2);
uc = mat2cell(u,siz);

% parallel step
if ~opt.inparallel
	for k=1:length(leaf)
		tic
		[z{k},J] = local_inverse(leaf{k},uc{k},in_border{k},diff{k},evalF,numcomp,Jac,tol_n/8);
		%        [l{k},u{k},p{k}] = lu(J,'vector');
		[L,U,p{k}] = lu(J,'vector');
		tim{k} = toc;
		fact{k} = tril(L,-1) + triu(U);
		z{k} = reshape(z{k},length(leaf{k}),opt.numcomp);
	end
else
	fprintf('  Starting parfor execution in SNK residual...')
	nc = opt.numcomp;
	parfor k=1:length(leaf)
		tic
		[z{k},J] = local_inverse(leaf{k},uc{k},in_border{k},diff{k},evalF,numcomp,Jac,tol_n/8);
		%        [l{k},u{k},p{k}] = lu(J,'vector');
		[L,U,p{k}] = lu(J,'vector');
		tim{k} = toc;
		fact{k} = tril(L,-1) + triu(U);
		z{k} = reshape(z{k},length(leaf{k}),nc);
	end
	fprintf('finished.\n')
end

%tim = cell2mat(tim)
z = cell2mat(z');
z = z(:);
luJ = struct('fact',fact,'p',p);

end


function [c,J] = local_inverse(approx,sol_k,border_k,diff_k,evalF,num_sols,Jac,tol_n)
% INPUT:
%      approx: leaf approximation
%       sol_k: given solution
%    border_k: border index for interface
%      diff_k: precomputed interface zone interpolation
%       evalF: residual function which returns Jacobian
%    num_sols: number of solutions
%  sol_length: length of solution
%
% OUTPUT
%           c: correction of solution
%          Jk: local Jocabian

% The residual is F(sol_k+z_k)
%            sol_k(border_k)+z_k(border_k)-B_k*u
%            (iterfacing at the zone interface)
	function F = residual(z)
		sol_length = length(approx);
		F = evalF(sol_k(:)+z,approx);
		F = reshape(F,sol_length,num_sols);
		z = reshape(z,sol_length,num_sols);
		F(border_k,:) = sol_k(border_k,:)+z(border_k,:) - diff_k;
		F = F(:);
	end

	function J = jac_fun(z)
		sol_length = length(approx);
		J = Jac(sol_k(:)+z(:),approx);
		E = eye(sol_length);
		
		for i=1:num_sols
			ind = false(sol_length*num_sols,1);
			ind((i-1)*sol_length+(1:sol_length)) = border_k;
			
			J(ind,:) = zeros(sum(ind),num_sols*sol_length);
			J(ind,(i-1)*sol_length+(1:sol_length)) = E(border_k,:);
		end
	end

	function [ sol,J ] = sol_and_jac( f,jac,u )
		sol = f(u);
		J = jac(u);
	end

init = zeros(numel(sol_k(:)),1);

options = optimoptions(@fsolve,'SpecifyObjectiveGradient',true,'MaxIterations',600,'FunctionTolerance',tol_n(1)/10,'Display','off');
[c,~,flag] = fsolve(@(u)sol_and_jac(@residual,@jac_fun,u),init,options);
if flag <=0
	warning('fsolve exited with code %i',flag)
end

J = jac_fun(c);

end




