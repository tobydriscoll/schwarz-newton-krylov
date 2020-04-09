function [ w ] = SNKJacobian(PUApprox,luJ,x,opt)
% SNKJacobian
% Matrix free evaluation for the Jacobian of the SNK method.
%
% INPUT:
%      PUApprox: PUApprox approximation
%         L,U,p: cell array of the LU decomposition of the local Jacobians.
%             x: solution
%
% NOTE x is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then x = [u1;u2;v1;v2].

numpoints = length(PUApprox);
leaf = PUApprox.leafArray;
numleaves = length(leaf);
x = reshape(x,numpoints,opt.numcomp);

% serial part
for k=1:numleaves
	degs{k} = leaf{k}.degs;
	[~,~,in_border,~] = FindBoundaryIndex2DSides(degs{k},PUApprox.leafArray{k}.domain,PUApprox.leafArray{k}.outerbox);
	
	z{k} = zeros(prod(degs{k}),opt.numcomp);
	
	if ~PUApprox.iscoarse
		z{k}(in_border,:) = leaf{k}.Binterp*x;
	else
		z{k}(in_border,:) = leaf{k}.CBinterp*x;
	end
	
	z{k} = z{k}(:);
end

fact = {luJ.fact};
p = {luJ.p};

xc = mat2cell(x,cellfun(@length,leaf));

% While parallelizable, hard to get real benefits here
if true %% ~opt.inparallel
	for k=1:length(PUApprox.leafArray)
		U = triu(fact{k});
		L = tril(fact{k},-1);  L(1:size(L,1)+1:end) = 1;
		v = U\(L\z{k}(p{k}));
		w{k} = reshape(v,length(leaf{k}),opt.numcomp)-xc{k};
	end
else
	fprintf('  Starting parfor execution in SNK Jacobian...')
	parfor k=1:length(PUApprox.leafArray)
		U = triu(fact{k});
		L = tril(fact{k},-1);  L(1:size(L,1)+1:end) = 1;
		v = U\(L\z{k}(p{k}));
		w{k} = reshape(v,length(leaf{k}),opt.numcomp)-xc{k};
	end
	fprintf('finished.\n')
end

w = cell2mat(w');
w = w(:);

