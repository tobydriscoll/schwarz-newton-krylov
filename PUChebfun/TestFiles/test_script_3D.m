u = [0.75 0.25 -0.75];

a = [25 25 25];

p = chebfunpref;
tol = p.cheb3Prefs.chebfun3eps;

 test_funs = {
             @(x,y,z) cos(2*pi*u(1)+a(1)*x+a(2)*y+a(3)*z);
             @(x,y,z)  1./(a(1)^-2+(x-u(1)).^2).*1./(a(2)^-2+(y-u(2)).^2).*1./(a(3)^-2+(z-u(3)).^2);
             @(x,y,z) 1./(1+a(1)*x+a(2)*y+a(3)*z).^(-3);
             @(x,y,z) exp(-a(1)*(x-u(1)).^2-a(2)*(y-u(2)).^2-a(3)*(z-u(3)).^2);
             @(x,y,z) 1./cosh(5*(x+y+z)).^2
             @(x,y,z) atan(5*(x+y)+z)};
        
CON_TIME_PU = zeros(length(test_funs),1);
INTERP_TIME_PU = zeros(length(test_funs),1);
INTERP_ERROR_PU = zeros(length(test_funs),1);
NUM_PTS_PU =  zeros(length(test_funs),1);

CON_TIME_CHEB  = zeros(length(test_funs),1);
INTERP_TIME_CHEB = zeros(length(test_funs),1);
INTERP_ERROR_CHEB = zeros(length(test_funs),1);
RANK_CHEB =  zeros(length(test_funs),1);

x = linspace(-1,1,100)';
G = {x x x};
[X,Y,Z] = ndgrid(x,x,x);

for i=1:length(test_funs)
    
    FV = test_funs{i}(X,Y,Z);
    
    M = max(abs(FV(:)));
    
    tic,TREE = PUchebfun(test_funs{i},[-1 1;-1 1;-1 1],'Degree',[129 129 129],'tol',tol); CON_TIME_PU(i) = toc;  
    tic,ef = TREE.evalfGrid(G);INTERP_TIME_PU(i) = toc;
    if any(isnan(ef(:)))
        error('found a nan');
    end
    E = abs(ef-FV); INTERP_ERROR_PU(i) = max(E(:))/M;
    NUM_PTS_PU(i) = length(TREE);
    
    tic, F = chebfun3(test_funs{i}); CON_TIME_CHEB(i) = toc;
    tic; ef = F(X,Y,Z); INTERP_TIME_CHEB(i) = toc;
    E = abs(ef-FV); INTERP_ERROR_CHEB(i) = max(E(:))/M;
    RANK_CHEB(i) = rank(F);
end

PU_TABLE = [INTERP_ERROR_PU CON_TIME_PU INTERP_TIME_PU NUM_PTS_PU];
CHEB_TABLE = [INTERP_ERROR_CHEB CON_TIME_CHEB INTERP_TIME_CHEB RANK_CHEB];
