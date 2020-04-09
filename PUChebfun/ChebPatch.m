classdef ChebPatch<LeafPatch
% ChebPatch PUFun class for representing n-d functions on lines, rectangles, and hyperectangles.
%
% This class represents a single tensor product polynomial The coefficients 
% are found through FFT.
 
% ChebPatch(var_struct) constructs a tensor product approximation 
% representing a function f, based on the options passed the structure
% var_struct.
%
% *var_struct.domain, the domain used for the Chebyshev polynomial.
%
% *var_struct.zone, the zone (non overlapping part from partition) used.
%
% *var_struct.outerbox the domain of the root of the tree.
%
% *var_struct.split_dim, an array of boolean indicies indicating if the 
%  approximation can be split in a given dimension.
%
% *var_struct.tol, the tolerance used for refinement: 'tol', 1e-b.
% 
% *var_struct.cdeg_in, the coarse degree to be used (if applicable) 
% 
% *var_struct.deg_in, the degree indices from the standard degrees in each 
% dimension for non square domains : 'degreeIndex', [ind_1,ind_2]. 

% Here the degrees can be chosen from the set [3 5 9 17 33 65 129].  
% So if 'degreeIndex', [5 5 5], the max degree of any approximate will be 
% 33 in each direction.
    properties
        swap_degs %temp holder for degs
        linOp
        linOp_f
        ClinOp
        Binterp
        CBinterp
        local_max
        orig_deg_in
        L
        U
        p
        split_all = false
    end
    
    properties (Access = protected)
        cdeg_in %index for the course degrees
        swap_deg_in
    end
    
    
    
    methods
        % Construct for the ChebPatch
        %
        %      Input:
        %       zone: (dim x 2) array indiciating array for the zone.
        %     domain: (dim x 2) array indiciating array for the domain.
        %   outerbox: (dim x 2) array indiciating array for the outerbox.
        %    degs_in: (dim x 1) integer array indicating the degree in each
        %                       dimension. Here the standard degrees are
        %                       [3 5 9 17 33 65 129].
        % split_flag: (dim x 1) boolean array indicating if the patch can
        %                       be split in any given dimension.
        %function obj = ChebPatch(domain,zone,outerbox,deg_in,split_flag,tol,cdeg_in)
        function obj = ChebPatch(var_struct)
            %Call superclass constructor
            obj = obj@LeafPatch(var_struct);
            obj.is_refined = false;
            obj.is_geometric_refined = true;
            
            if isfield(var_struct, 'split_all')
                obj.split_all = var_struct.split_all;
            end
            
        end
            
        %Returns structure of parameters
        function p_struct = params(obj)
            p_struct.outerbox = obj.outerbox;
            p_struct.zone = obj.zone;
            p_struct.domain = obj.domain;
            p_struct.split_flag = obj.split_flag;
            p_struct.tol = obj.tol;
            p_struct.cdegs = obj.cdegs;
            p_struct.degs = obj.degs;
            p_struct.orig_degs = obj.orig_degs;
            p_struct.orig_cdegs = obj.orig_cdegs;
            p_struct.split_all = obj.split_all;
        end
        
        % Sets the values to be used for interpolation
        %
        % Input:
        %     f: values sampled at obj.points.
        function [Max] = sample(obj,f,grid_opt,~)
            
            if(nargin==2)
                grid_opt = false;
            end
            
            %Just assume we sample f for right now.
            if ~isnumeric(f)
                
                if ~grid_opt
                    if obj.dim==1
                        obj.coeffs = chebtech2.vals2coeffs(f(obj.points()));
                    elseif obj.dim==2
                        points = obj.points();
                        V = reshape(f(points(:,1),points(:,2)),obj.degs);
                        obj.coeffs = chebfun2.vals2coeffs(V);
                    else
                        points = obj.points();
                        V = reshape(f(points(:,1),points(:,2),points(:,3)),obj.degs);
                        obj.coeffs = chebfun3.vals2coeffs(V); 
                    end
                else
                    %If a function is more efficient on a grid, evaluate it
                    %as so.
                    V = f(obj.leafGrids());
                    if obj.dim==2
                        obj.coeffs = chebfun2.vals2coeffs(V);
                    else
                        obj.coeffs = chebfun3.vals2coeffs(V);
                    end
                end
            else
                switch obj.dim
                    case 1
                        [~,n2] = size(f);
                        if n2 == 1
                            V = f';
                        else
                            V = f;
                        end
                        obj.coeffs = chebtech2.vals2coeffs(V);
                    case 2
                        V = reshape(f,obj.degs);
                        obj.coeffs = chebfun2.vals2coeffs(V);
                    otherwise
                        V = reshape(f,obj.degs);
                        obj.coeffs = chebfun3.vals2coeffs(V);
                end
            end
            
            obj.values = V;
            
            Max = max(abs(V(:)));
            obj.local_max = Max;
            
        end
        
        % The method determines if a splitting is needed, and creates
        % the new children if splitting is required.
        %
        %     Input:
        %   overlap: overlap intended to be used for the splitting
        %
        %    Output:
        %     Child: obj if no splitting is needed. If a splitting
        %            is needed, Child is the PUPatch object with
        %            the new children.
        function Child = splitleaf(obj,Max,set_vals)
            
            if nargin==2
                set_vals = false;
            end
            
            vscale = Max;
            
            loc_tol = obj.tol^(7/8);
            
            cutoff = zeros(obj.dim,1);
            isHappy = false(obj.dim,1);
            
            
            if obj.dim==1
                fCol = chebtech2({[],obj.coeffs});
                hscale = diff(obj.domain);
                tol = loc_tol*max(vscale./obj.local_max,hscale);
                cutoff(1) = length(simplify(fCol, tol))+1;
                
                p = chebfunpref;
                p.chebfuneps = obj.tol*vscale;
                
                [isHappy(1), cutoff(1)] = happinessCheck(fCol, [], [], [], p);
            else
                
                for k=1:obj.dim
                    
                    if obj.split_flag(k)

                        colChebtech = chebfun3t.unfold(obj.coeffs, k);
                        colChebtech = sum(abs(colChebtech),2);
                        fCol = chebtech2({[],colChebtech});
                        hscale = diff(obj.domain(k,:));
                        
                        %p = chebfunpref;
                        %p.chebfuneps = loc_tol;
                        
                        tol = loc_tol*max(vscale./obj.local_max,hscale);
                        
                        %[isHappy(k), cutoff(k)] = happinessCheck(fCol, [], [], [], p);
                        
                        cutoff(k) = length(simplify(fCol, tol))+1;
                        
                        isHappy(k) = cutoff(k)<obj.degs(k);
                        
                    end
                    
                end
            end
            
            for k=1:obj.dim
                if isHappy(k)
                    obj.degs(k) = cutoff(k);
                    if obj.degs(k)<obj.cdegs(k)
                        obj.cdegs(k) = obj.degs(k);
                        obj.values = [];
                    end
                    obj.split_flag(k) = false;
                end
            end
            
            if obj.dim==1
                obj.coeffs = obj.coeffs(1:obj.degs);
            elseif obj.dim==2
                obj.coeffs = obj.coeffs(1:obj.degs(1),1:obj.degs(2));
            else
                obj.coeffs = obj.coeffs(1:obj.degs(1),1:obj.degs(2),1:obj.degs(3));
            end
            
            if ~any(obj.split_flag)
                %The leaf is refined, so return it.
                obj.is_refined = true;
                obj.cheb_length = prod(obj.degs);
                Child = obj;
                if set_vals
                    if obj.dim==1
                        obj.values = chebtech2.coeffs2vals(obj.coeffs);
                    elseif obj.dim==2
                        obj.values = chebfun2.coeffs2vals(obj.coeffs);
                    else
                        obj.values = chebfun3.coeffs2vals(obj.coeffs);
                    end
                end
            else
                
                
                Child = obj;
                %Go through and split in each unresolved direction
                for k=1:obj.dim
                    if obj.split_flag(k) || obj.split_all
                        if Child.is_leaf
                            Child = split(obj,k,set_vals);
                        else
                            Child.split(k,set_vals);
                        end
                    end
                end
                
                
                
            end
        end
        
        
        % The method determines will split a child into along
        % a dimension.
        %
        %     Input:
        %   overlap: overlap intended to be used for the splitting
        %
        %    Output:
        %     Child: the PUPatch with the two new children.
        function Child = split(obj,split_dim,set_vals)
            
            if nargin == 2
                set_vals = false;
            end
            
            Children = cell(1,2);
            
            %The width of the overlap
            delta = 0.5*obj.overlap*diff(obj.zone(split_dim,:));
            
            zone0 = obj.zone;
            zone1 = obj.zone;
            
            region0 = obj.domain;
            region1 = obj.domain;
            
            m = sum(obj.zone(split_dim,:))/2;
            
            zone0(split_dim,:) = [obj.zone(split_dim,1) m];
            zone1(split_dim,:) = [m obj.zone(split_dim,2)];
            
            region0(split_dim,:) = [max(obj.outerbox(split_dim,1),obj.zone(split_dim,1)-delta) m+delta];
            region1(split_dim,:) = [m-delta,min(obj.outerbox(split_dim,2),obj.zone(split_dim,2)+delta)];
            
            
            struct0 = obj.params;
            struct0.domain = region0; struct0.zone = zone0;
            
            struct1 = obj.params;
            struct1.domain = region1; struct1.zone = zone1;
            
            
            Children{1} = ChebPatch(struct0);
            Children{2} = ChebPatch(struct1);
            
            Children{1}.orig_degs = obj.orig_degs;
            Children{1}.orig_deg_in = obj.orig_deg_in;
            
            Children{2}.orig_degs = obj.orig_degs;
            Children{2}.orig_deg_in = obj.orig_deg_in;
            
            pu_domain = obj.domain;
            
            pu_domain(split_dim,:) = [region0(split_dim,1) region1(split_dim,2)];
            
            %Return the PUPatch with the new children
            Child = PUPatch(pu_domain,obj.zone,Children,split_dim);
            
            if set_vals
                for k=1:2
                    V = obj.evalfGrid(Child.children{k}.leafGrids());
                    Child.children{k}.sample(V);
                    Child.children{k}.values = V;
                end
            end
        end
        
        % String function for object.
        %
        %     Input:
        %
        %     Child: string of object.
        function str = toString(obj)
            str = '';
            for i=1:obj.dim
                str = strcat(str,sprintf(' [%0.3f %0.3f]',obj.domain(i,1),obj.domain(i,2)));
                if i<obj.dim
                    str = strcat(str,' x ');
                end
            end
            str = strcat(str,' degrees: ');
            for i=1:obj.dim
                str = strcat(str,sprintf(' %d ',obj.degs(i)));
            end
            
            str = strcat(str,sprintf(' length %d', obj.length));
        end
        
        function IsGeometricallyRefined = IsLeafGeometricallyRefined(obj)
            IsGeometricallyRefined = true;
        end
        
        function V = Getvalues(obj)
            V = obj.values(:);
        end
        
        
        function reset(obj)
            obj.is_refined = false;
            obj.is_geometric_refined = false;
            obj.degs = obj.orig_degs;
            obj.coeffs = zeros(obj.degs);
            obj.cdegs = obj.orig_cdegs;
            obj.values = zeros(obj.degs);
            obj.split_flag = ones(obj.dim,1);
        end
       
        function Child = splitleafGeom(obj)
            Child = obj;
        end
    end
end