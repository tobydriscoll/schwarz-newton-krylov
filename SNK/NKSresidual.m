function [z] = NKSresidual(sol,PUApprox,f,opt)
% PDE residual, including boundary/interface differences.
%
% INPUT
%     sol: given solution at patches
%
%     PUApprox: PUApprox approximation
%
%     f: residual function f(x,p) for solution x and local
%        approximation p. f(x,p) evaluates the residual on the domain
%        of p.
%
% OUTPUT
%     z: residual of solution / identity at inner boundary of patches

% NOTE sol is presumed to be ordered by solution components first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two components u1 v1, u2 v2. Then sol = [u1;u2;v1;v2].

sol = reshape(sol,length(PUApprox),opt.numcomp);

% Figure out offset index for each patch
leaf = PUApprox.leafArray;
step = zeros(length(leaf),1);
for k=2:length(leaf)
	step(k) = step(k-1) + length(leaf{k-1});
end

for k=1:length(leaf)	
	degs = leaf{k}.degs;
	
	% Get logical indices for interface and boundary points. Output is
	% given for all indices, as well as the indices along each of the
	% sides.
	[~,~,interface{k}] = FindBorders(degs,leaf{k}.domain,leaf{k}.outerbox);
	
	index{k} = sol(step(k)+(1:prod(degs)),:);
	
	if ~PUApprox.iscoarse
		diff{k} = leaf{k}.Binterp*sol;
	else
		diff{k} = leaf{k}.CBinterp*sol;
	end
end

% parallel step
if ~opt.inparallel
	for k=1:length(leaf)
		% Assume z is of the form [u1 u2 ... un]
		z{k} = local_residual(leaf{k},index{k},interface{k}.all,diff{k},f,opt.numcomp);
	end
else
	fprintf('  Starting parallel execution in NKS residual...')
	nc = opt.numcomp;
	parfor k=1:length(leaf)
		% Assume z is of the form [u1 u2 ... un]
		z{k} = local_residual(leaf{k},index{k},interface{k}.all,diff{k},f,nc);
	end
	fprintf('finished.\n')
end

z = cell2mat(z');
z = z(:);

end

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

function F = local_residual(approx,sol_k,border_k,diff_k,evalF,num_sols)

F = evalF(sol_k(:),approx);
sol_length = length(approx);
F = reshape(F,sol_length,num_sols);
F(border_k,:) = sol_k(border_k,:) - diff_k;

end

