classdef (Abstract) PUfun < handle & matlab.mixin.Copyable
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
    % *The domain of the function: 'domain' , [a_1 b_1;a_2 b_2;...;a_d b_d]
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
        ChebRoot
        leafArray
    end
    
    
    methods (Abstract)
        addTree = plus(obj,Tree2);
        subTree = minus(obj,Tree2);
        MultTree = mtimes(obj,Tree2);
        DivTree = mrdivide(obj,Tree2);
        PowTree = mpower(obj,p);
    end
    
    methods
        
        %This method will split leaves if they are determined to be
        %unrefined.
        function splitleaves(obj,Max,set_vals)
            if obj.ChebRoot.is_leaf
                obj.ChebRoot = obj.ChebRoot.splitleaf(Max,set_vals);
            else
                obj.ChebRoot.PUsplit(Max,set_vals);
            end
            obj.clean();
        end
        
        
        % ef = evalf(obj,X)
        % This method evalutes the PUM approximation at X.
        %Input:
        %   x          : this constructs matrices given f|x.
        %Output:
        %   ef         : vector of approximation at X
        function ef = evalf(obj,X)
            ef = obj.ChebRoot.evalf(X);
        end
        
        % ef = evalf(obj,X)
        % This method evalutes the PUM approximation at at grid X.
        %Input:
        %   x          : this constructs matrices given f|x.
        %Output:
        %   ef         : grid of approximation at X
        function ef = evalfGrid(obj,X)
            ef = obj.ChebRoot.evalfGrid(X);
        end
        
        function Max = sample(obj,f)
            Max = obj.ChebRoot.sample(f);
        end
        
        function clean(obj)
            obj.ChebRoot.clean();
            obj.leafArray = obj.ChebRoot.collectLeaves();
        end
        
        % int = sum(obj)
        % This method computes the integral of obj
        function int = sum(obj)
            
            int = 0;
            for i=1:length(obj.leafArray)
                
                if ~isempty(obj.leafArray{i}.zone)
                    X = cell(1,obj.ChebRoot.dim);
                    
                    for j=1:obj.ChebRoot.dim
                        [X{j},W{j}] = chebpts(obj.leafArray{i}.degs(j),obj.leafArray{i}.zone(j,:));
                    end
                    
                    vals = obj.leafArray{i}.evalfGrid(X);
                    
                    if obj.ChebRoot.dim==1
                        int = int + W{1}*vals;
                    elseif obj.ChebRoot.dim==2
                        int = int + W{2}*(W{1}*vals).';
                    else
                        int = int + W{2}*(W{1}*chebfun3.txm(vals,W{3},3))';
                    end
                    
                end
            end
        end
        
        function N = norm(obj)
            N = sum(obj^2);
        end
        
        % diff_Tree = diff(obj,diff_dim,order)
        % This method computes thd approximation of the derivative
        %Input:
        %Output:
        %   diff_Tree  : PU approximation of derivative
        function diff_Tree = diff(obj,diff_dim,order)

            diff_Tree = copy(obj);
            
            for i=1:length(obj.leafArray)
                diff_Tree.leafArray{i}.coeffs = DiffCoeffs(obj.leafArray{i},diff_dim,order);
            end
            
        end
        
        % div_Tree = div(obj)
        % This method computes thd approximation of the derivative
        %Input:
        %Output:
        %   div_Tree  : PU approximation of divergence
        function div_Tree = div(obj)
            div_Tree = copy(obj);
            
            for i=1:length(obj.leafArray)
                if obj.ChebRoot.dim==1
                        div_Tree.leafArray{i}.coeffs = DiffCoeffs(obj.leafArray{i},1,1);
                elseif obj.ChebRoot.dim==2
                        div_Tree.leafArray{i}.coeffs = DiffCoeffs(obj.leafArray{i},1,1)+DiffCoeffs(obj.leafArray{i},2,1);
                else
                        div_Tree.leafArray{i}.coeffs = DiffCoeffs(obj.leafArray{i},1,1)+DiffCoeffs(obj.leafArray{i},2,1)+DiffCoeffs(obj.leafArray{i},3,1);
                end
            end
            
        end
        
        % lap_Tree = lap(obj)
        % This method computes thd approximation of the Laplacian
        %Input:
        %Output:
        %   lap_Tree  : PU approximation of Laplacian
        function lap_Tree = lap(obj)
            lap_Tree = copy(obj);
            
            for i=1:length(obj.leafArray)
                if obj.ChebRoot.dim==1
                        lap_Tree.leafArray{i}.coeffs = DiffCoeffs(obj.leafArray{i},1,2);
                elseif obj.ChebRoot.dim==2
                        lap_Tree.leafArray{i}.coeffs = DiffCoeffs(obj.leafArray{i},1,2)+DiffCoeffs(obj.leafArray{i},2,2);
                else
                        lap_Tree.leafArray{i}.coeffs = DiffCoeffs(obj.leafArray{i},1,2)+DiffCoeffs(obj.leafArray{i},2,2)+DiffCoeffs(obj.leafArray{i},3,2);
                end
            end
        end
        
        % ln = length(obj)
        % This returns the length of the tree
        function ln = length(obj)
            ln = length(obj.ChebRoot);
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
    
    methods (Static)
        %This method slices an array along a certain dimension
        function out = slice(A, ix, dim)
            subses = repmat({':'}, [1 ndims(A)]);
            subses{dim} = ix;
            out = A(subses{:});
        end
        
        
        % T_add = add(T_1,T_2,add_f)
        % This method adds T_1 and T_2
        %Input:
        %   T_1,T_2    : trees to be added together
        %   add_f      : function for adding T_1 and T_2 using PU
        %Output:
        %   T_add      : PU approximation of the sum
        function T_merge = fast_add(T_1,T_2,T_merge,split_dim)
            
            if T_1.is_leaf && T_2.is_leaf
                T_merge.degs = max(T_1.degs,T_2.degs);
                T_merge.orig_degs = T_merge.degs;
                    T_merge.coeffs = zeros(T_merge.degs);
                    
                    if isequal(T_1.domain,T_2.domain)
                        if T_merge.dim==1
                            T_merge.coeffs(1:T_1.degs) = T_1.coeffs;
                            T_merge.coeffs(1:T_2.degs) = T_merge.coeffs(1:T_2.deg)+ T_2.coeffs;
                        elseif T_merge.dim==2
                            T_merge.coeffs(1:T_1.degs(1),1:T_1.degs(2)) = T_1.coeffs;
                            T_merge.coeffs(1:T_2.degs(1),1:T_2.degs(2)) = T_merge.coeffs(1:T_2.degs(1),1:T_2.degs(2))+ T_2.coeffs;
                        elseif T_merge.dim==3
                            T_merge.coeffs(1:T_1.degs(1),1:T_1.degs(2),1:T_1.degs(3)) = T_1.coeffs;
                            T_merge.coeffs(1:T_2.degs(1),1:T_2.degs(2),1:T_2.degs(3)) = T_merge.coeffs(1:T_2.degs(1),1:T_2.degs(2),1:T_2.degs(3))+ T_2.coeffs;
                        end
                    else
                        
                        V = T_1.evalfGrid(T_merge.leafGrids())+T_2.evalfGrid(T_merge.leafGrids());
                        
                        if T_merge.dim==1
                            T_merge.coeffs = chebtech2.vals2coeffs(V);
                        elseif T_merge.dim==2
                            T_merge.coeffs = chebfun2.vals2coeffs(V);
                        elseif T_merge.dim==3
                            T_merge.coeffs = chebfun3.vals2coeffs(V);
                        end
                        
                    end
                
            elseif T_1.is_leaf && ~T_2.is_leaf
                T_merge = T_merge.split(T_2.splitting_dim);
                
                for k=1:2
                    if ~T_merge.children{k}.is_null
                        T_merge.children{k} = PUchebfun.fast_add(T_1,T_2.children{k},T_merge.children{k},T_merge.splitting_dim);
                    end
                end
                
            elseif ~T_1.is_leaf && T_2.is_leaf
                T_merge = T_merge.split(T_1.splitting_dim);
                
                for k=1:2
                    if ~T_merge.children{k}.is_null
                        T_merge.children{k} = PUchebfun.fast_add(T_1.children{k},T_2,T_merge.children{k},T_merge.splitting_dim);
                    end
                end
                
            elseif T_1.splitting_dim == T_2.splitting_dim
                T_merge = T_merge.split(T_1.splitting_dim);
                
                for k=1:2
                    if ~T_merge.children{k}.is_null
                        T_merge.children{k} = PUchebfun.fast_add(T_1.children{k},T_2.children{k},T_merge.children{k},T_merge.splitting_dim);
                    end
                end
                
            else
                
                if split_dim==0
                    next_split=1;
                else
                    next_split = mod(split_dim,T_1.dim)+1;
                end
                
                while true
                    if T_1.splitting_dim == next_split
                        First_T = T_1;
                        Last_T = T_2;
                        break;
                    elseif T_2.splitting_dim == next_split
                        First_T = T_2;
                        Last_T = T_1;
                        break;
                    else
                        next_split = mod(next_split,T_1.dim)+1;
                    end
                end
                
                
                if ~Last_T.split_flag(First_T.splitting_dim)
                    T_merge = T_merge.split(First_T.splitting_dim);
                    
                    for k=1:2
                        if ~T_merge.children{k}.is_null
                            T_merge.children{k} = PUchebfun.fast_add(First_T.children{k},Last_T,T_merge.children{k},T_merge.splitting_dim);
                        end
                    end
                    
                else
                    T_merge = T_merge.refine(@(x)T_1.evalfGrid(x)+T_2.evalfGrid(x),true);
                end
                
                
            end
        end
        
        % T_add = add(T_1,T_2,add_f)
        % This method adds T_1 and T_2
        %Input:
        %   T_1,T_2    : trees to be added together
        %   add_f      : function for adding T_1 and T_2 using PU
        %Output:
        %   T_add      : PU approximation of the sum
        function T_merge = fast_subtract(T_1,T_2,T_merge,split_dim)
            
            if T_1.is_leaf && T_2.is_leaf
                T_merge.degs = max(T_1.degs,T_2.degs);
                T_merge.orig_degs = T_merge.degs;
                T_merge.coeffs = zeros(T_merge.degs);
                    
                    if isequal(T_1.domain,T_2.domain)
                        if T_merge.dim==1
                            T_merge.coeffs(1:T_1.degs) = T_1.coeffs;
                            T_merge.coeffs(1:T_2.degs) = T_merge.coeffs(1:T_2.deg) - T_2.coeffs;
                        elseif T_merge.dim==2
                            T_merge.coeffs(1:T_1.degs(1),1:T_1.degs(2)) = T_1.coeffs;
                            T_merge.coeffs(1:T_2.degs(1),1:T_2.degs(2)) = T_merge.coeffs(1:T_2.degs(1),1:T_2.degs(2)) - T_2.coeffs;
                        elseif T_merge.dim==3
                            T_merge.coeffs(1:T_1.degs(1),1:T_1.degs(2),1:T_1.degs(3)) = T_1.coeffs;
                            T_merge.coeffs(1:T_2.degs(1),1:T_2.degs(2),1:T_2.degs(3)) = T_merge.coeffs(1:T_2.degs(1),1:T_2.degs(2),1:T_2.degs(3)) - T_2.coeffs;
                        end
                    else
                        
                        V = T_1.evalfGrid(T_merge.leafGrids())-T_2.evalfGrid(T_merge.leafGrids());
                        
                        if T_merge.dim==1
                            T_merge.coeffs = tech.vals2coeffs(V);
                        elseif T_merge.dim==2
                            T_merge.coeffs = chebfun2.vals2coeffs(V);
                        elseif T_mergre.dim==3
                            T_merge.coeffs = chebfun3.vals2coeffs(V);
                        end
                        
                    end
                
            elseif T_1.is_leaf && ~T_2.is_leaf
                T_merge = T_merge.split(T_2.splitting_dim);
                
                for k=1:2
                    if ~T_merge.children{k}.is_null
                        T_merge.children{k} = PUchebfun.fast_subtract(T_1,T_2.children{k},T_merge.children{k},T_merge.splitting_dim);
                    end
                end
                
            elseif ~T_1.is_leaf && T_2.is_leaf
                T_merge = T_merge.split(T_1.splitting_dim);
                
                for k=1:2
                    if ~T_merge.children{k}.is_null
                        T_merge.children{k} = PUchebfun.fast_subtract(T_1.children{k},T_2,T_merge.children{k},T_merge.splitting_dim);
                    end
                end
                
            elseif T_1.splitting_dim == T_2.splitting_dim
                T_merge = T_merge.split(T_1.splitting_dim);
                
                for k=1:2
                    if ~T_merge.children{k}.is_null
                        T_merge.children{k} = PUchebfun.fast_subtract(T_1.children{k},T_2.children{k},T_merge.children{k},T_merge.splitting_dim);
                    end
                end
                
            else
                
                if split_dim==0
                    next_split=1;
                else
                    next_split = mod(split_dim,T_1.dim)+1;
                end
                
                while true
                    if T_1.splitting_dim == next_split
                        First_T = T_1;
                        Last_T = T_2;
                        break;
                    elseif T_2.splitting_dim == next_split
                        First_T = T_2;
                        Last_T = T_1;
                        break;
                    else
                        next_split = mod(next_split,T_1.dim)+1;
                    end
                end
                
                
                if ~Last_T.split_flag(First_T.splitting_dim)
                    T_merge = T_merge.split(First_T.splitting_dim);
                    
                    for k=1:2
                        if ~T_merge.children{k}.is_null
                            T_merge.children{k} = PUchebfun.fast_subtract(First_T.children{k},Last_T,T_merge.children{k},T_merge.splitting_dim);
                        end
                    end
                    
                else
                    T_merge = T_merge.refine(@(x)T_1.evalfGrid(x)-T_2.evalfGrid(x),true,true);
                end
                
                
            end
        end
        
        % T_add = add(T_1,T_2,add_f)
        % This method adds T_1 and T_2
        %Input:
        %   T_1,T_2    : trees to be added together
        %   add_f      : function for adding T_1 and T_2 using PU
        %Output:
        %   T_add      : PU approximation of the sum
        function T_merge = fast_multiply(T_1,T_2,T_merge,split_dim)
            
            if T_1.is_leaf && T_2.is_leaf
                
                T_merge.degs = max(T_1.degs,T_2.degs);
                T_merge.orig_degs = T_merge.degs;
                T_merge = T_merge.refine(@(x)T_1.evalfGrid(x).*T_2.evalfGrid(x),true,true);
                
            elseif T_1.is_leaf && ~T_2.is_leaf
                T_merge = T_merge.split(T_2.splitting_dim);
                
                for k=1:2
                    if ~T_merge.children{k}.is_null
                        T_merge.children{k} = PUchebfun.fast_multiply(T_1,T_2.children{k},T_merge.children{k},T_merge.splitting_dim);
                    end
                end
                
            elseif ~T_1.is_leaf && T_2.is_leaf
                T_merge = T_merge.split(T_1.splitting_dim);
                
                for k=1:2
                    if ~T_merge.children{k}.is_null
                        T_merge.children{k} = PUchebfun.fast_multiply(T_1.children{k},T_2,T_merge.children{k},T_merge.splitting_dim);
                    end
                end
                
            elseif T_1.splitting_dim == T_2.splitting_dim
                T_merge = T_merge.split(T_1.splitting_dim);
                
                for k=1:2
                    if ~T_merge.children{k}.is_null
                        T_merge.children{k} = PUchebfun.fast_multiply(T_1.children{k},T_2.children{k},T_merge.children{k},T_merge.splitting_dim);
                    end
                end
                
            else
                
                if split_dim==0
                    next_split=1;
                else
                    next_split = mod(split_dim,T_1.dim)+1;
                end
                
                while true
                    if T_1.splitting_dim == next_split
                        First_T = T_1;
                        Last_T = T_2;
                        break;
                    elseif T_2.splitting_dim == next_split
                        First_T = T_2;
                        Last_T = T_1;
                        break;
                    else
                        next_split = mod(next_split,T_1.dim)+1;
                    end
                end
                
                
                if ~Last_T.split_flag(First_T.splitting_dim)
                    T_merge = T_merge.split(First_T.splitting_dim);
                    
                    for k=1:2
                        if ~T_merge.children{k}.is_null
                            T_merge.children{k} = PUchebfun.fast_multiply(First_T.children{k},Last_T,T_merge.children{k},T_merge.splitting_dim);
                        end
                    end
                    
                else
                    T_merge = T_merge.refine(@(x)T_1.evalfGrid(x).*T_2.evalfGrid(x),true,true);
                end
                
            end
        end
        
        % T_add = add(T_1,T_2,add_f)
        % This method adds T_1 and T_2
        %Input:
        %   T_1,T_2    : trees to be added together
        %   add_f      : function for adding T_1 and T_2 using PU
        %Output:
        %   T_add      : PU approximation of the sum
        function T_merge = fast_divide(T_1,T_2,T_merge,split_dim)
            
            if T_1.is_leaf && T_2.is_leaf
                
                T_merge.degs = max(T_1.degs,T_2.degs);
                T_merge.orig_degs = T_merge.degs;
                T_merge = T_merge.refine(@(x)T_1.evalfGrid(x)./T_2.evalfGrid(x),true,true);
                
            elseif T_1.is_leaf && ~T_2.is_leaf
                T_merge = T_merge.split(T_2.splitting_dim);
                
                for k=1:2
                    if ~T_merge.children{k}.is_null
                        T_merge.children{k} = PUchebfun.fast_divide(T_1,T_2.children{k},T_merge.children{k},T_merge.splitting_dim);
                    end
                end
                
            elseif ~T_1.is_leaf && T_2.is_leaf
                T_merge = T_merge.split(T_1.splitting_dim);
                
                for k=1:2
                    if ~T_merge.children{k}.is_null
                        T_merge.children{k} = PUchebfun.fast_divide(T_1.children{k},T_2,T_merge.children{k},T_merge.splitting_dim);
                    end
                end
                
            elseif T_1.splitting_dim == T_2.splitting_dim
                T_merge = T_merge.split(T_1.splitting_dim);
                
                for k=1:2
                    if ~T_merge.children{k}.is_null
                        T_merge.children{k} = PUchebfun.fast_divide(T_1.children{k},T_2.children{k},T_merge.children{k},T_merge.splitting_dim);
                    end
                end
                
            else
                
                if split_dim==0
                    next_split=1;
                else
                    next_split = mod(split_dim,T_1.dim)+1;
                end
                
                while true
                    if T_1.splitting_dim == next_split
                        First_T = T_1;
                        Last_T = T_2;
                        break;
                    elseif T_2.splitting_dim == next_split
                        First_T = T_2;
                        Last_T = T_1;
                        break;
                    else
                        next_split = mod(next_split,T_1.dim)+1;
                    end
                end
                
                
                if ~Last_T.split_flag(First_T.splitting_dim)
                    T_merge = T_merge.split(First_T.splitting_dim);
                    
                    for k=1:2
                        if ~T_merge.children{k}.is_null
                            T_merge.children{k} = PUchebfun.fast_divide(First_T.children{k},Last_T,T_merge.children{k},T_merge.splitting_dim);
                        end
                    end
                    
                else
                    T_merge = T_merge.refine(@(x)T_1.evalfGrid(x)./T_2.evalfGrid(x),true,true);
                end
                
            end
        end
        
        
        % T_power = power(Tree,p)
        % This method computes Tree^p
        %Input:
        %   Tree       : trees to be powered
        %   p          : power to use
        %Output:
        %   T_power    : PU approximation of the power
        function T_power = power(Tree,p)
            T_power = copy(Tree);
            
            T_power.reset();
            
            T_Array = Tree.collectLeaves();
            leafArray = T_power.collectLeaves();
            
            
            for i=1:length(T_Array)
                leafArray{i} = refine(leafArray{i},@(x)T_Array{i}.evalfGrid(x).^p,true,true);
            end
        end
    end
end

