classdef (Abstract) Patch < handle & matlab.mixin.Copyable
% Patch PUFun abstract class for representing a patch. Here a patch is
% defined to be either a PU approximation or a single tensor product
% polynomial.

    properties
        domain
        zone
        boundaryIndex
        cheb_length
        Root
        is_leaf
        is_refined = false
        is_geometric_refined = false
        iscoarse = false;
        dim
        tol
        outerbox
        split_flag %array indicating if we will split along a dimension
        is_null = false;
    end
    
    properties (Constant)
        overlap = 0.08;
    end
    
    methods (Abstract)
        points = points(obj)
        ef = evalf(obj,x,dim,order)
        ef = evalfGrid(obj,x,dim,order)
        ln = length(obj)
        sample(obj,f)
        plotdomain(obj)
    end
    
    methods
        
        %Returns logical indicies of points in the domain 
        function domain_ind = InDomain(obj,x)
            [num_pts,~] = size(x);
            domain_ind = true(num_pts,1);
            for i=1:obj.dim
                domain_ind = domain_ind & ...
                    ( x(:,i)>= obj.domain(i,1) & x(:,i)<= obj.domain(i,2));
            end
        end
        
        %Returns logical indicies of points in the zone 
        function domain_ind = InZone(obj,x)
            domain_ind = true(size(x,1),1);
            for i=1:obj.dim
                domain_ind = domain_ind & ...
                    ( x(:,i)>= obj.zone(i,1) & x(:,i)<= obj.zone(i,2));
            end
        end
        

        %determine which points along the grid belong to the domain
        function [sub_grid,split_ind] = IndDomainGrid(obj,x)
            sub_grid = cell(1,obj.dim);
            split_ind = cell(1,obj.dim);
            
            for i=1:obj.dim
                ind = x{i}>=obj.domain(i,1) & x{i}<=obj.domain(i,2);
                
                split_ind{i} = ind;
                
                sub_grid{i} = x{i}(ind);
            end
        end
        
        %2D method for showing patches
        function show(obj,level,step)
            % Make a pretty graph showing the domains (2D)
            assert(obj.dim==2 || obj.dim==1,'Must be 1-D or 2-D')
            if nargin==1  % user call
                level = 0;
                step = 1;
            elseif nargin==2
                step = 1;
            end
            
            LEAVES = obj.collectLeaves();
            
            for i=1:length(LEAVES)
                if obj.dim==1
                    hold on;
                    plot(LEAVES{i}.domain,(level)*ones(1,2),'LineWidth',2,'color','black');
                    hold off;
                elseif obj.dim==2
                    x = LEAVES{i}.domain(1,[1 1 2 2]);
                    y = LEAVES{i}.domain(2,[1 2 2 1]);
                    z = level*ones(1,4);
                    patch(x,y,z,rand(1,3),'facealpha',0.2,'edgecolor','none')
                end
                level = level+step;
            end
            
        end
        
        %This method will recursively determine the length, as well as
        %determine the splitting structure
        function clean(obj)
            if obj.is_leaf
                %obj.split_flag = false(size(obj.split_flag));
                obj.cheb_length = prod(obj.degs);
            else
                
                clean(obj.children{1});
                clean(obj.children{2});
                obj.split_flag = obj.children{1}.split_flag | obj.children{2}.split_flag;
                obj.split_flag(obj.splitting_dim) = true;
                
                for i=1:obj.dim
                    obj.domain(i,:) = [min(obj.children{1}.domain(i,1),min(obj.children{2}.domain(i,1))) max(obj.children{1}.domain(i,2),min(obj.children{2}.domain(i,2)))];
                end
                
                obj.cheb_length = obj.children{1}.cheb_length+obj.children{2}.cheb_length;
            end
        end
        
        %  collectLeaves(obj,leaves)
        %  Recursive function that collects leaves into cell array.
        %
        %     Input:
        %    leaves: current list of leaves
        %
        %    Output:
        %    LEAVES: leaves list with children of patch added.
        function [LEAVES]= collectLeaves(obj)
            if obj.is_leaf
                LEAVES = {obj};
            else
                LEAVES = collectLeaves_recurse(obj,{});
            end
        end
        
        %this function recursively collects the patches
        function [LEAVES]= collectLeaves_recurse(obj,leaves)
            for k=1:2
                if obj.children{k}.is_leaf && ~obj.children{k}.is_null
                    leaves{length(leaves)+1} = obj.children{k};
                elseif ~obj.children{k}.is_null
                    leaves = obj.children{k}.collectLeaves_recurse(leaves);
                end
            end
            LEAVES = leaves;
        end
        
        %This method, given a function, will refine the patch.
        function fun_obj = refine(obj,f,grid_opt,fast_opt)
            
            if nargin<3
                grid_opt = false;
                fast_opt = false;
            elseif nargin<4
                fast_opt = false;
            end
            
            fun_obj = obj;
            
            while ~fun_obj.is_refined
                
                Max = fun_obj.sample(f,grid_opt,fast_opt);
                
                if fun_obj.is_leaf
                    fun_obj = obj.splitleaf(Max);
                else
                    fun_obj.PUsplit(Max);
                end
                
            end
            
        end
   
        
    end
end

