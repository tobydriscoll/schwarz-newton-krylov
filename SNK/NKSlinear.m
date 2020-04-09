function [z] = NKSlinear(sol,PUApprox,J,opt)
% Residual of the linearized PDE.
%
% INPUT:
%      sol: given solution at patches
%      PUApprox: PUApprox approximation
%      J: Jcell array of the local Jacobians.
%
% OUTPUT:
%          z: residual of solution, identity at inner boundary of patches
%
% NOTE sol is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then sol = [u1;u2;v1;v2].

leaf = PUApprox.leafArray;
numleaves = length(leaf);
step = zeros(numleaves,1);
sol = reshape(sol,length(PUApprox),opt.numcomp);

% Figure out starting index for each patch
for k=2:numleaves
	step(k) = step(k-1) + length(leaf{k-1});
end

for k=1:numleaves
	degs = leaf{k}.degs;
	
	%This function returns the logical indicies of the gamma and outer
	%boundry interface. Out put is given for all indicies, as well as the
	%indicies along each of the sides
	[~,~,in_border{k},~] = FindBoundaryIndex2DSides(degs,leaf{k}.domain,leaf{k}.outerbox);

	sol_loc{k} = sol(step(k)+(1:prod(degs)),:);
	
	if ~PUApprox.iscoarse
		diffr{k} = leaf{k}.Binterp*sol;
	else
		diffr{k} = leaf{k}.CBinterp*sol;
	end
end


% parallel step
if ~opt.inparallel
	for k=1:numleaves
		%Assume z is of the form [u1 u2 ... un]
		[z{k}] = local_residual(leaf{k},sol_loc{k},in_border{k},diffr{k},J{k},opt.numcomp);
	end
else
	fprintf('  Starting parallel execution in NKS linear...')
	nc = opt.numcomp;
	parfor k=1:numleaves
		%Assume z is of the form [u1 u2 ... un]
		[z{k}] = local_residual(leaf{k},sol_loc{k},in_border{k},diffr{k},J{k},nc);
	end
	fprintf('finished.\n')
end

z = cell2mat(z');
z = z(:);

end

function F = local_residual(approx,sol_k,border_k,diff_k,J,nc)
% INPUT:
%      approx: leaf approximation
%       sol_k: given solution
%    border_k: border index for interface
%      diff_k: precomputed interface zone interpolation
%       evalF: residual function which returns Jacobian
%    opt.numcomp: number of solutions
%  sol_length: length of solution
%
% OUTPUT
%           c: correction of solution
%          Jk: local Jocabian

F = J*sol_k(:);
sol_length = length(approx);
F = reshape(F,sol_length,nc);
F(border_k,:) = sol_k(border_k,:) - diff_k;

end

