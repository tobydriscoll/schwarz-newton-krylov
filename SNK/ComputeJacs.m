function [J,l,u,p] = ComputeJacs(sol,PUApprox,Jac)
% Computes Jacobian for NKS method
%
% INPUT:
%      PUApprox: PUApprox approximation
%         sol: given solution
%       evalF: residual function which returns Jacobian
%    num_sols: number of solutions
%
% OUTPUT:
%          z: correction of solution
%          J: cell array of local Jacobians
%
% NOTE sol is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then sol = [u1;u2;v1;v2].

num_sols = length(sol)/length(PUApprox);
step = zeros(length(PUApprox.leafArray),1);
sol = reshape(sol,length(PUApprox),num_sols);

%Figure out starting index for each patch
for k=2:length(PUApprox.leafArray)
	step(k) = step(k-1) + length(PUApprox.leafArray{k-1});
end

for k=1:length(PUApprox.leafArray)
	
	degs = PUApprox.leafArray{k}.degs;
	
	%This will be sol_length*num_sols
	sol_loc{k} = sol(step(k)+(1:prod(degs)),:);
	
	%This function returns the logical indicies of the gamma and outer
	%boundry interface. Out put is given for all indicies, as well as the
	%indicies along each of the sides
	[~,~,in_border{k},~] = FindBoundaryIndex2DSides(degs,PUApprox.leafArray{k}.domain,PUApprox.leafArray{k}.outerbox);
	
end

%parallel step

leafs = PUApprox.leafArray;

for k=1:length(leafs)
	[J{k},l{k},u{k},p{k}] = local_Jac(leafs{k},sol_loc{k},in_border{k},num_sols,Jac);
end

end

function [J,l,u,p] = local_Jac(approx,sol_k,border_k,num_sols,Jac)
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

J = jac_fun(sol_k(:));
[l,u,p] = lu(J,'vector');

	function J = jac_fun(z)
		sol_length = length(approx);
		J = Jac(z(:),approx);
		E = eye(sol_length);
		
		for i=1:num_sols
			ind = false(sol_length*num_sols,1);
			ind((i-1)*sol_length+(1:sol_length)) = border_k;
			
			J(ind,:) = zeros(sum(ind),num_sols*sol_length);
			J(ind,(i-1)*sol_length+(1:sol_length)) = E(border_k,:);
		end
	end

end





