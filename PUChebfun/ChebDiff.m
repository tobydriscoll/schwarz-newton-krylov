function [ Dcx ] = ChebDiff(N,p,dom)
% This function computes a differentiation matrices for a n-1 degree
% Chebyshev polynomial.

Dcx =zeros(n,n);

for i=1:floor(n/2)
    Dcx(n+1+2*(i-1)*n:n+1:end) = (2+4*(i-1)):2:(n-1)*2;
end

Dcx(1,:) = 0.5*Dcx(1,:);

end

