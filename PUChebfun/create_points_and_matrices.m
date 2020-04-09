N=3;

M=7;

chebpoints = cell(M,1);

chebmatrices = cell(M,2);

chebweights = cell(M,1);

for i=1:M
    chebpoints{i} = chebpts(N);
    chebmatrices{i,1} = diffmat(N,1);
    chebmatrices{i,2} = diffmat(N,2);
    chebweights{i} = chebtech2.barywts(N);
    N = N+(N-1);
end