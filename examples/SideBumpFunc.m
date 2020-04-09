function [ f ] = SideBumpFunc(x,domain,m)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

c = domain(2)-m;

f = zeros(size(x));

ind = x<domain(2) & x>c;
x(ind) = (x(ind)-domain(2))/m;
f(ind) = exp(-(x(ind).^2)./(1-x(ind).^2));

f(x>=domain(2)) = 1;

