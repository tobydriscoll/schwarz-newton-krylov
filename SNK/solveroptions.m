function opt = solveroptions(varargin)

ip = inputParser;
ip.KeepUnmatched = true;
ip.addParameter('reltol',1e-5,@isnumeric);
ip.addParameter('abstol',1e-6,@isnumeric);
ip.addParameter('inparallel',false,@islogical);
ip.addParameter('griddepth',0,@isnumeric);
ip.addParameter('coarsetol',1e-4);
ip.addParameter('twolevel',true,@islogical);
ip.addParameter('numwork',2,@isnumeric);
ip.addParameter('numcomp',1,@isnumeric);
ip.addParameter('solver',@SNKsolver);
parse(ip,varargin{:});

opt = ip.Results;
end
