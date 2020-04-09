u = [0.75 0.25];

a = [5 10];

[~,F1] = cheb.gallery2('squarepeg');
[~,F2] = cheb.gallery2('tiltedpeg');

p = chebfunpref;
tol = p.cheb2Prefs.chebfun2eps;

test_funs = {@(x,y) log(1+10^5*(x.^2+y.^2));
            @(x,y) atan(10^2*(x+y.^2));
            @(x,y) 10^(-4)./((10^(-4)+x.^2).*(10^(-4)+y.^2));
            @(x,y) franke(x,y);
            @(x,y) cos(2*pi*u(1)+a(1)*x+a(2)*y);
            @(x,y)  1./(a(1)^-2+(x-u(1)).^2).*1./(a(2)^-2+(y-u(2)).^2);
            @(x,y) 1./(1+a(1)*x+a(2)*y).^(-3);
            @(x,y) exp(-a(1)*(x-u(1)).^2-a(2)*(y-u(2)).^2);
            F1;
            F2};
        
        
CON_TIME_PU = zeros(length(test_funs),1);
INTERP_TIME_PU = zeros(length(test_funs),1);
INTERP_ERROR_PU = zeros(length(test_funs),1);
NUM_PTS_PU =  zeros(length(test_funs),1);

CON_TIME_CHEB  = zeros(length(test_funs),1);
INTERP_TIME_CHEB = zeros(length(test_funs),1);
INTERP_ERROR_CHEB = zeros(length(test_funs),1);
RANK_CHEB =  zeros(length(test_funs),1);

x = linspace(-1,1,100)';
G = {x x};
[X,Y] = ndgrid(x,x);

for i=1:length(test_funs)
    tic,TREE = PUchebfun(test_funs{i},[-1 1;-1 1],'Degree',[129 129],'tol',tol); CON_TIME_PU(i) = toc;
    tic,ef = TREE.evalfGrid(G);INTERP_TIME_PU(i) = toc;
    
    if any(isnan(ef(:)))
        error('found a nan');
    end
    
    
    FV = test_funs{i}(X,Y);
    M = max(abs(FV(:)));
    
    E = abs(ef-FV); INTERP_ERROR_PU(i) = max(E(:))/M;
    NUM_PTS_PU(i) = length(TREE);
    
    tic, F = chebfun2(test_funs{i}); CON_TIME_CHEB(i) = toc;
    tic; ef = F(X,Y); INTERP_TIME_CHEB(i) = toc;
    E = abs(ef-FV); INTERP_ERROR_CHEB(i) = max(E(:))/M;
    RANK_CHEB(i) = rank(F);
end

PU_TABLE = [INTERP_ERROR_PU CON_TIME_PU INTERP_TIME_PU NUM_PTS_PU];
CHEB_TABLE = [INTERP_ERROR_CHEB CON_TIME_CHEB INTERP_TIME_CHEB RANK_CHEB];
