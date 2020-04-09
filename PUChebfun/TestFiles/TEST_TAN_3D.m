clear;

c = 0.2;

NT = 20;
NP = 20;


AVT = 1;

degs = [6 6 6];

domain = [-1 1;-1 1;-1 1];

ang_t = linspace(0,pi/2,NT);
ang_p = linspace(0,pi/2,NT);

k = 10;

f2t = @(x,y,z,t,p) atan(k*(sin(p)*cos(t)*x+sin(p)*sin(t)*y+cos(p)*z));
ft = @(x,t,p) f2t(x(:,1),x(:,2),x(:,3),t,p);

TIMES1 = zeros(NT,NP);
TIMES2 = zeros(NT,NP);

for i=1:NT
    for k=1:NP
    AV = zeros(AVT,1);
    
    for j=1:AVT
        tic;
        TREE = PUchebfun(@(x,y,z)f2t(x,y,z,ang_t(i),ang_p(k)),'domain',domain,'degreeIndex',degs);
        AV(j)=toc;
    end
    
    TIMES1(i,k) = mean(AV);
    
    for j=1:AVT
        tic;
        F = chebfun3(@(x,y,z)f2t(x,y,z,ang_t(i),ang_p(k)));
        AV(j)=toc;
    end
    
    TIMES2(i,k) = mean(AV);
    a=1;
    end
end

