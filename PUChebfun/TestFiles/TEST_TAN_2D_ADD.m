clear;

c = 0.2;

NT = 10;

AVT = 1;

degs = [7 7];

domain = [-1 1;-1 1];

ang = linspace(0,pi/2,NT);

k = 2.5e2;

g2t = @(x,y) atan(k*(x+y));
gt = @(x) g2t(x(:,1),x(:,2));

f2t = @(x,y,t) atan(k*(cos(t)*x+sin(t)*y));
ft = @(x,t) f2t(x(:,1),x(:,2),t);

TIMES1 = zeros(NT,1);
TIMES2 = zeros(NT,1);

TREE1 = PUFun(domain,degs,gt,1e-12);
F1 = chebfun2(g2t);

for i=1:NT
    
    AV = zeros(AVT,1);
    
    TREE2 = PUFun(domain,degs,@(x)ft(x,ang(i)),1e-12);
    
    for j=1:AVT
        tic;
        TREE3 = TREE1*TREE2;
        AV(j)=toc;
    end
    
    TIMES1(i) = mean(AV);
    
    F2 = chebfun2(@(x,y)f2t(x,y,ang(i)));
    
    for j=1:AVT
        tic;
        F = F1.*F2;
        AV(j)=toc;
    end
    
    TIMES2(i) = mean(AV);
    
    a=1;
    
end

