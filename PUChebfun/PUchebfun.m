classdef PUchebfun < PUfun
% PUchebfun   PUchebfun class for representing n-d functions on lines, rectangles
% and cuboids using partition of unitities.
% 
% This class represents smooth multivariate functions on hypercubes up to 
% dimension 3 with a partition of unity approximation. This class 
% automatically finds a set of overlapping domains that are adapted to the 
% features of the function, and the blends locally defined Chebyshev 
% approximations on the domains with a partition of unity.
% 
% PUchebfun(f) constructs a partition of unity approximation representing f on 
% the domain [-1 1]^k, where k is the number of variables of f 
% (i.e. the dimension of the domain of f). Functions must be vectorized.
% 
% PUchebfun(f,[a_1 b_1;a_2 b_2;...;a_d b_d]) constructs a partition of unity 
% approximation representing f on the domain [a_1 b_1;a_2 b_2;...;a_d b_d]. 
% Functions must be vectorized. 
% 
% PUchebfun(f,[a_1 b_1;a_2 b_2;...;a_d b_d],varargin) constructs a partition of unity approximation 
% representing f, based on the options passed into with varargin; that is 
% PUchebfun(f,'perf1',perf1,'pref2',pref2,..) is called. This preferences that 
% can be set are:
% 
% 
% *The degree indices from the standard degrees in each dimension : 
% 'degreeIndex', [ind_1,ind_2, ... ind_d]. 
% 
% % *The tolerance used for adaptation : 
% 'tol', tol. 
%
% Here the degrees can be chosen from the set [3 5 9 17 33 65 129].  
% So if 'degreeIndex', [5 5 5], the max degree of any approximate will be 
% 33 in each direction. 

    properties
        degs
        tol
        domain = [];
        grid_opt = false;
        iscoarse = false;
    end
    
    methods
        
        %function obj = PUchebfun(domain,deg_in,f,tol,grid_opt)
        function obj = PUchebfun(varargin)
            
            is_Patch = false;
                
                if length(varargin)==1
                    
                    if isa(varargin{1},'Patch')
                        obj.ChebRoot = varargin{1};
                        is_Patch = true;
                        obj.leafArray = obj.ChebRoot.collectLeaves();
                        obj.domain = obj.ChebRoot.domain;
                    else
                        f = varargin{1};
                        obj.domain = repmat([-1 1],nargin(f),1);
                        cheb_struct.domain = obj.domain;
                        obj.ChebRoot = ChebPatch(cheb_struct);
                    end
                    

                    
                elseif length(varargin)==2
                    f = varargin{1};
                    obj.domain = varargin{2};
                    cheb_struct.domain = obj.domain;
                    obj.ChebRoot = ChebPatch(cheb_struct);
                    
                else
                    f = varargin{1};
                    obj.domain = varargin{2};
                    cheb_struct.domain = obj.domain;
                    varargin(1:2) = [];
                    args = varargin;
                    while ( ~isempty(args) )
                        if strcmpi(args{1}, 'gridOption')
                            obj.grid_opt = args{2};
                        elseif strcmpi(args{1}, 'Degree')
                            obj.degs = args{2};
                            cheb_struct.orig_degs = args{2};
                            cheb_struct.degs = args{2};
                        elseif strcmpi(args{1}, 'SplitAll')
                            cheb_struct.split_all = args{2};    
                        elseif strcmpi(args{1}, 'tol')
                            obj.tol = args{2};
                            cheb_struct.tol = args{2};
                        elseif strcmpi(args{1}, 'CoarseDegree')
                            cheb_struct.orig_cdegs = args{2};
                            cheb_struct.cdegs = args{2};
                        else
                            error(strcat(args{1},' is not a valid parameter.'));
                        end
                        args(1:2) = [];
                    end
                    
                    obj.ChebRoot = ChebPatch(cheb_struct);
                    
                end
                
                
                obj.tol = obj.ChebRoot.tol;
                
                
                if ~is_Patch
                    refine(obj,f);
                end
    end
            
        
        % refine(obj,f,grid_opt)
        % This method refines the tree by f(x).
        %Input:
        %   f          : the function to split on
        %   grid_opt   : boolean value indicating if
        %                function is evaluated for grids;
        %                must take cell array of grids
        function refine(obj,f)
            
            
            while ~obj.ChebRoot.is_refined
                
                Max = obj.ChebRoot.sample(f,obj.grid_opt);
                
                if obj.ChebRoot.is_leaf
                    obj.ChebRoot = obj.ChebRoot.splitleaf(Max);
                else
                    obj.ChebRoot.PUsplit(Max);
                end
                
            end
            
            
            if ~obj.ChebRoot.is_leaf
                obj.ChebRoot.findIndex([]);
                obj.leafArray = obj.ChebRoot.collectLeaves();
            else
                obj.leafArray = {obj.ChebRoot};
            end
            
            obj.ChebRoot.clean();
            
        end
        
        
        function Patch = newRoot(obj)
            vars.domain = obj.domain;
            vars.degs = obj.degs;
            Patch = ChebPatch(vars);
        end
        
        function v = Getvalues(obj)
            v = obj.ChebRoot.Getvalues(); 
        end
        
        
        % addTree = plus(obj,Tree2)
        % This method adds obj and Tree2
        %Input:
        %   Tree2      : the other tree
        %Output:
        %   addTree    : new tree of the sum
        function addTree = plus(obj,Tree2)
            
            add_T = newRoot(obj);

            addTreeRoot = PUfun.fast_add(obj.ChebRoot,Tree2.ChebRoot,add_T,0);
            
            addTreeRoot.clean();
            
            addTree = PUchebfun(addTreeRoot);
            
        end
        
        % subTree = minus(obj,Tree2)
        % This method subtracts obj and Tree2
        %Input:
        %   Tree2      : the other tree
        %Output:
        %   subTree    : new tree of the difference
        function subTree = minus(obj,Tree2)
            
            sub_T = newRoot(obj);
            
            subTreeRoot = PUfun.fast_subtract(obj.ChebRoot,Tree2.ChebRoot,sub_T,0);
            
            subTreeRoot.clean();
            
            subTree = PUchebfun(subTreeRoot);
        end
        
        % Coarsen(obj)
        % This method Coarsens each of the patches
        function reset(obj)
            obj.ChebRoot.reset();
            obj.ChebRoot.clean();
        end
        
        % MultTree = mtimes(obj,Tree2)
        % This method multiplies obj and Tree2
        %Input:
        %   Tree2      : the other tree
        %Output:
        %   MultTree   : new tree of the product
        function MultTree = mtimes(obj,Tree2)
            
            mult_T = newRoot(obj);
            
            multTreeRoot = PUfun.fast_multiply(obj.ChebRoot,Tree2.ChebRoot,mult_T,0);
            
            multTreeRoot.clean();
            
            MultTree = PUchebfun(multTreeRoot);
            
        end
        
        % DivTree = mrdivide(obj,Tree2)
        % This method divides obj and Tree2
        %Input:
        %   Tree2      : the other tree
        %Output:
        %   DivTree    : new tree of the quotient
        function DivTree = mrdivide(obj,Tree2)
            
            div_T = newRoot(obj);
            
            DivTreeRoot = PUfun.fast_divide(obj.ChebRoot,Tree2.ChebRoot,div_T,0);
            
            DivTreeRoot.clean();
            
            DivTree = PUchebfun(DivTreeRoot);
            
            
        end
        
        % PowTree = mpower(obj,p)
        % This method computes obj to the power p
        %Input:
        %   p          : the power to be used
        %Output:
        %   PowTree    : the new tree of the power
        function PowTree = mpower(obj,p)
            
            PowTreeRoot = PUfun.power(obj.ChebRoot,p);
            
            PowTree = PUchebfun(PowTreeRoot);
            
        end
        
        % Coarsen(obj)
        % This method Coarsens each of the patches
        function Coarsen(obj)
            obj.ChebRoot.Coarsen();
            obj.iscoarse = true;
            obj.clean();
        end
        
        function cvals = Fine2Coarse(obj,vals,k)
            
            if 2==nargin
                k=0;
            end
            
            cvals = cell(length(obj.leafArray),1);

            step = zeros(length(obj.leafArray),1);
            
            %Figure out starting index for each patch
            for j=2:length(obj.leafArray)
                step(j) = step(j-1) + length(obj.leafArray{j-1});
            end


            for i=1:length(obj.leafArray)
                sub_ind = (1:length(obj.leafArray{i}))+step(i);
                cvals{i} = obj.leafArray{i}.Fine2Coarse(vals(sub_ind),k);
            end
            
            cvals = cell2mat(cvals);
            
        end
        
        function rvals = Coarse2Fine(obj,vals)
            
           % rvals = obj.ChebRoot.Coarse2Fine(vals);
            
                        
            rvals = cell(length(obj.leafArray),1);

            step = zeros(length(obj.leafArray),1);
            
            %Figure out starting index for each patch
            for j=2:length(obj.leafArray)
                step(j) = step(j-1) + length(obj.leafArray{j-1});
            end


            for i=1:length(obj.leafArray)
                sub_ind = (1:length(obj.leafArray{i}))+step(i);
                rvals{i} = obj.leafArray{i}.Coarse2Fine(vals(sub_ind));
            end
            
            rvals = cell2mat(rvals);
            
        end
        
        % Refines(obj)
        % This method Refines each of the patches
        function Refine(obj)
            obj.ChebRoot.Refine();
            obj.iscoarse = false;
            obj.clean();
        end
        
        
        % ln = length(obj)
        % This returns the length of the tree
        function ln = length(obj)
            ln = length(obj.ChebRoot);
        end
        
        % disp(obj)
        % Returns string of object
        function disp(obj)
            disp(obj.ChebRoot.toString());
        end
        
        function Setvalues(obj,f)
            obj.ChebRoot.Setvalues(f);
        end
        
            
        % disp(obj)
        % Returns color plot of patches
        function show(obj,level,step)
            if nargin==1  % user call
                level = 0;
                step = 1;
            elseif nargin==2
                step = 1;
            end
            
            show(obj.ChebRoot,level,step)
        end
        
        % plot(obj)
        % Returns a plot of the function in 2D
        function plot(obj)
            
            assert(obj.ChebRoot.dim ==2 || obj.ChebRoot.dim ==1,'Must be 1-D or 2-D');
            
            if obj.ChebRoot.dim ==1
                domain = obj.ChebRoot.domain;
                x = linspace(domain(1,1),domain(1,2),100)';
                ef = obj.evalf(x);
                plot(x,ef,'LineWidth',2,'color','blue');
            else
                domain = obj.ChebRoot.domain;
                x = linspace(domain(1,1),domain(1,2),100)';
                y = linspace(domain(2,1),domain(2,2),100)';
                ef = obj.evalfGrid({x y});
                [X,Y] = ndgrid(x,y);
                defaultOpts = {'facecolor', 'flat', 'edgealpha', .5, 'edgecolor', 'none'};
                surf(X,Y,ef,defaultOpts{:});
                xlabel('x');
                ylabel('y');
            end
        end
    end
    
    methods(Access = protected)
        % Override copyElement method:
        function cpObj = copyElement(obj)
            % Make a shallow copy of all properties
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            
            cpObj.ChebRoot = obj.ChebRoot.copy();
            
            cpObj.leafArray = cpObj.ChebRoot.collectLeaves();
        end
    end
     
        
end

