if opt.inparallel
    pp = gcp('nocreate');
    if isempty(pp) || pp.NumWorkers ~=opt.numwork
        delete(pp);
        pp = parpool(opt.numwork);
    end

    bstart = pp.ticBytes;
    tstart = tic;
    [ sol,stats ] = opt.solver(f,Jac,init,F,opt);
    stats.elapsedtime = toc(tstart);
    stats.tocbytes = pp.tocBytes(bstart);
else
    tstart = tic;
    [ sol,stats ] = opt.solver(f,Jac,init,F,opt);
    stats.elapsedtime = toc(tstart);
end

opt
stats

% Transfer last solution component to tree structure
solr = reshape(sol,length(F),opt.numcomp);
F.sample(solr(:,1));
%plot(F);  % uncomment to plot 