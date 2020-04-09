classdef PUPatch<Patch
% LSPatch2D PUFun class for representing a patch with a single splitting.

% This class is used as part of a recursively defined tree, where
% objects at the leaves are represented with the LeafPatch.m class.
%
% LSPatch2D(domain,zone,children,splitting_dim) will construct the
% PUPatch. Here the domain is the overlaping domain (tyically the bounding
% box of the domains of the children), zone is the nonoverlapping domain
% from the partition, children is an cell array of the two children, and
% splitting_dim is the dimension the patch is split in.
    properties
        children
        splitting_dim
        index = [];
    end
    
    properties (Constant)
        weights = PUWeights();
    end
    
    methods
        % Construct for the PUPatch
        %
        %        Input:
        %        domain: (dim x 2) array indiciating array for the domain.
        %          zone: (dim x 2) array indiciating array for the zone.
        %      children: (2 x 1) cell array of Patch objects;
        % splitting_dim: dimension the patch is split in.
        function obj = PUPatch(domain,zone,children,splitting_dim)
            
            if isempty(children{1}) && ismepty(children{2})
                error('one child must not be empty');
            end
            
            obj.outerbox = children{1}.outerbox;
            obj.zone = zone;
            obj.domain = domain;
            [obj.dim,~] = size(obj.domain);
            obj.cheb_length = children{1}.length()+children{2}.length();
            obj.children = children;
            obj.splitting_dim = splitting_dim;
            obj.is_leaf = false;
            obj.is_refined = false;
            obj.split_flag = obj.children{1}.split_flag | obj.children{2}.split_flag;
        end
        
        % Returns the length of the patch
        %   Output:
        %       ln: total number of interpolating points
        function ln = length(obj)
            ln = obj.cheb_length;
        end
        
        % Returns cell array of grids on leaves
        function grids = leafGrids(obj)
            gridsl = leafGrids(obj.children{1});
            gridsr = leafGrids(obj.children{2});
            if obj.children{1}.is_leaf && obj.children{2}.is_leaf
                grids = {gridsl gridsr};
            elseif obj.children{1}.is_leaf && ~ obj.children{2}.is_leaf
                grids = {gridsl gridsr{:}};
            elseif obj.children{2}.is_leaf
                grids = {gridsl{:} gridsr};
            else
                grids = {gridsl{:} gridsr{:}};
            end
        end
        
        % Returns vector of interpolating points for the leaf.
        function pts = points(obj)
            pts = [obj.children{1}.points();obj.children{2}.points()];
        end
        
        % Will split the children of the patch if they are unrefined.
        %
        %     Output:
        % is_refined: returns if the patch is refined.
        function is_refined = PUsplit(obj,Max,set_vals)
            
            if nargin ==2
                set_vals = false;
            end
                
            is_refined = true;
            is_geom_refined = true;
            
            
            for k=1:2
                if ~obj.children{k}.is_null && ~obj.children{k}.is_refined
                    if obj.children{k}.is_leaf
                        obj.children{k} = obj.children{k}.splitleaf(Max,set_vals);
                        is_refined = obj.children{k}.is_refined && is_refined;
                        is_geom_refined = obj.children{k}.is_geometric_refined && is_geom_refined;
                    else
                        obj.children{k}.PUsplit(Max,set_vals);
                        is_refined = is_refined && obj.children{k}.is_refined;
                        is_geom_refined = is_geom_refined && obj.children{k}.is_geometric_refined;
                    end
                end
            end
            
            obj.is_geometric_refined = is_geom_refined;
            obj.is_refined = is_refined;
        end
        
                % Will split the children of the patch if they are unrefined.
        %
        %     Output:
        % is_refined: returns if the patch is refined.
        function is_geometric_refined = PUsplitGeom(obj)

            is_geometric_refined = true;
            
            for k=1:2
                if obj.children{k}.is_leaf && ~obj.children{k}.is_null
                    obj.children{k} = obj.children{k}.splitleafGeom();
                    is_geometric_refined = obj.children{k}.is_geometric_refined && is_geometric_refined;
                else
                    obj.children{k}.PUsplitGeom();
                    is_geometric_refined =  obj.children{k}.is_geometric_refined && is_geometric_refined;
                end
            end
            
            obj.is_geometric_refined = is_geometric_refined;
            
        end
        
        % This method samples the leaves of the patch.
        %
        %     Input:
        %         f: vector of values for the interpolating values on the
        %            leaves depth first, or an anonymous function.
        function [Max] = sample(obj,f,grid_opt,fast_opt)
            
            if isnumeric(f) || ~obj.is_refined
                
                if nargin==2
                    grid_opt = false;
                    fast_opt = false;
                elseif nargin==3
                    fast_opt = false;
                end
                    
                
                MaxC = [-inf -inf];
                
                if ~isnumeric(f)
                    for k=1:2
                        if ~obj.children{k}.is_null
                            MaxC(k) = obj.children{k}.sample(f,grid_opt,fast_opt);
                        end
                    end
                else
                    
                    if obj.children{1}.is_null
                        MaxC(2) = obj.children{2}.sample(f,grid_opt,fast_opt);
                    elseif obj.children{2}.is_null
                        MaxC(1) = obj.children{1}.sample(f,grid_opt,fast_opt);
                    else
                    end
                    MaxC(1) = obj.children{1}.sample(f(1:length(obj.children{1})),grid_opt,fast_opt);
                    MaxC(2) = obj.children{2}.sample(f(length(obj.children{1})+1:end),grid_opt,fast_opt);
                end
                Max = max(MaxC);
                
            else
                
                Max = -inf;
                
            end
            
        end
        
        
        %  findIndex(obj,index)
        %  This function finds the index (i.e. the path from the root to
        %  the patch recursively for all the leaves of this patch.
        %
        %Input:
        %   index    : the index from the root to this patch
        %Output:
        function findIndex(obj,index)
            
            for k=1:2
                index_c = [index k];
                if obj.children{k}.is_leaf
                    obj.children{k}.index = index_c;
                else
                    obj.children{k}.findIndex(index_c);
                end
            end
            
        end
        
        %  evalfGrid(obj,X)
        %  Evaluates the approximation on a grid.
        %
        %Input:
        %    X: cell array of grid values.
        function vals = evalfGrid(obj,X)
            
            if obj.dim==1
                vals = evalf(obj,X{1});
            else
                [sum,dotprod] = obj.evalfGrid_recurse(X);
                vals = dotprod./sum;
            end
        end
        
        %  evalf(obj,X)
        %  Evaluates the approximation on a list of points.
        %
        %Input:
        %    X: list of points.
        function vals = evalf(obj,X)
            [sum,dotprod] = obj.evalf_recurse(X);
            vals = dotprod./sum;
        end
        
        %  evalf_recurse(obj,X)
        %  Evaluates the values need recursively for the approximation on a
        %  list of points X.
        %
        %Input:
        %    X: list of points.
        function [sum,dotprod] = evalf_recurse(obj,X)
            
            [num_pts,~] = size(X);
            
            dotprod = zeros(num_pts,1);
            
            sum = zeros(num_pts,1);
            
            
            %calculate values for the children
            for k=1:2
                ind = obj.children{k}.InDomain(X);
                
                
                if any(ind) && ~obj.children{k}.is_null
                    
                    if ~obj.children{k}.is_leaf
                        [sumk,dotprodk] = obj.children{k}.evalf_recurse(X(ind,:));
                    else
                        sumk = obj.children{k}.evalfBump(X(ind,:));
                        dotprodk = sumk.*obj.children{k}.evalf(X(ind,:));
                    end
                    
                    sum(ind) = sum(ind) + sumk;
                    dotprod(ind) = dotprod(ind) + dotprodk;
                end
                
            end
        end
        
        %  evalfGrid_recurse(obj,X)
        %  Evaluates the values need recursively for the approximation on a
        %  grid X.
        %
        %Input:
        %    X: cell array of grid values.
        function [sum,dotprod] = evalfGrid_recurse(obj,X)
            
            grid_lengths = cellfun(@(x)size(x,1),X);
            
            sum = zeros(grid_lengths);
            
            dotprod = zeros(grid_lengths);
            
            for k=1:2
                [sub_grid{k},sub_ind{k}] = obj.children{k}.IndDomainGrid(X);
            end
            
            for i=1:obj.dim
                com_ind{i} = sub_ind{1}{i} & sub_ind{2}{i};
                
                for k=1:2
                    com_sub_ind{k}{i} = com_ind{i}(sub_ind{k}{i});
                end
            end
            
            
            %            calculate values for the children
            
            for k=1:2
                if all(cellfun(@any,sub_ind{k})) && ~obj.children{k}.is_null
                    if ~obj.children{k}.is_leaf
                        [sumk,dotprodk] = obj.children{k}.evalfGrid_recurse(sub_grid{k});
                    else
                        sumk = obj.children{k}.evalfGridBump(sub_grid{k});
                        dotprodk = sumk.*obj.children{k}.evalfGrid(sub_grid{k});
                    end
                    
                    if obj.dim == 2
                        
                        sum(sub_ind{k}{1},sub_ind{k}{2}) = sum(sub_ind{k}{1},sub_ind{k}{2})+sumk;
                        dotprod(sub_ind{k}{1},sub_ind{k}{2}) = dotprod(sub_ind{k}{1},sub_ind{k}{2})+dotprodk;
                        
                    else
                        
                        if k==1
                            sum(sub_ind{k}{1},sub_ind{k}{2},sub_ind{k}{3}) = sumk;
                            dotprod(sub_ind{k}{1},sub_ind{k}{2},sub_ind{k}{3}) = dotprodk;
                            
                        else
                            
                        sum(sub_ind{k}{1},sub_ind{k}{2},sub_ind{k}{3}) = sum(sub_ind{k}{1},sub_ind{k}{2},sub_ind{k}{3})+sumk;
                        dotprod(sub_ind{k}{1},sub_ind{k}{2},sub_ind{k}{3}) = dotprod(sub_ind{k}{1},sub_ind{k}{2},sub_ind{k}{3})+dotprodk;
                        
                        end
                        
                    end
                end
            end
        end
        
        %  interpSparseMatrixZone(obj,X)
        %  Creates a sparse interpolating matrix for a list of points X.
        %
        % Output:
        %      M: sparse interpolating matrix.
        function M = interpSparseMatrixZone(obj,X)
            [ii,jj,zz] = obj.interpMatrixZone_vecs(X);
            M = sparse(ii,jj,zz,size(X,1),length(obj));
        end
        
        %  interpSparseMatrix(obj,X)
        %  Creates a sparse interpolating matrix for a list of points X.
        %
        % Output:
        %      M: sparse interpolating matrix.
        function M = interpSparseMatrix(obj,X)
            [ii,jj,zz,sum] = obj.interpMatrix_vecs(X);
            zz = zz./sum(ii);
            M = sparse(ii,jj,zz,size(X,1),length(obj));
        end
        
        
        %  interpSparseMatrixZoneGrid(obj,X)
        %  Creates a sparse interpolating matrix along the zones
        %  for a grid X.
        %
        % Output:
        %      M: sparse interpolating matrix.
        function M = interpSparseMatrixZoneGrid(obj,X)
            [ii,jj,zz] = obj.interpMatrixZoneGrid_vecs(X);
            grid_lengths = cellfun(@(x)length(x),X);
            M = sparse(ii,jj,zz,prod(grid_lengths),length(obj));
        end
        
        
        %  interpSparseMatrixZoneGrid(obj,X)
        %  Creates a sparse interpolating matrix along the zones
        %  for a grid X.
        %
        % Output:
        %      M: sparse interpolating matrix.
        function M = interpSparseMatrixGrid(obj,X)
            [ii,jj,zz,sum] = obj.interpMatrixGrid_vecs(X);
            sum = sum(:);
            zz = zz./sum(ii);
            grid_lengths = cellfun(@(x)length(x),X);
            M = sparse(ii,jj,zz,prod(grid_lengths),length(obj));
        end
        
        
        %  interpSparseMatrixZone_vecs(obj,X)
        %  Creates a sparse interpolating matrix along the zones
        %  for a list of points X.
        %
        %     Input:
        %      grid: cell array of grid values.
        %
        %    Output:
        %  ii,jj,kk: vectors used to construct the sparse matrix.
        function [ii,jj,zz] = interpMatrixZoneGrid_vecs(obj,grid)
            
            ii = []; jj = []; zz = [];
            grid_lengths = cellfun(@(x)length(x),grid);
            
            [n,~] = size(grid_lengths);
            
            if n>1
                grid_lengths = grid_lengths';
            end
            
            sub_grids = cell(2,1);
            
            inds = cell(2,1);
            
            not_empty_child = true(2,1);
            
            sub_grids{1} = grid;
            sub_grids{2} = grid;
            
            midpoint = mean(obj.zone(obj.splitting_dim,:));
            
            inds{1} = grid{obj.splitting_dim} < midpoint;
            inds{2} = ~inds{1};
            
            for k=1:2
                sub_grids{k}{obj.splitting_dim} = sub_grids{k}{obj.splitting_dim}(inds{k});
                not_empty_child(k) = all(cellfun(@(x)~isempty(x),sub_grids{k}));
            end
            
            in_index = (1:prod(grid_lengths))';
            
            step = 0;
            
            for k=1:2
                if not_empty_child(k)
                    if obj.children{k}.is_leaf && ~obj.children{k}.is_null
                        [ii_k,jj_k,zz_k] = find(obj.children{k}.interpMatrixGrid(sub_grids{k}));
                        
                        if size(ii_k,2)~=1
                            ii_k = ii_k';
                        end
                        
                        if size(jj_k,2)~=1
                            jj_k = jj_k';
                        end
                        
                        if size(zz_k,2)~=1
                            zz_k = zz_k';
                        end
                        
                    else
                        [ii_k,jj_k,zz_k] = obj.children{k}.interpMatrixZoneGrid_vecs(sub_grids{k});
                    end
                    
                    if obj.dim==2
                        ind_k = false(grid_lengths(1),grid_lengths(2));
                        
                        if obj.splitting_dim==1
                            ind_k(inds{k},:) = true;
                        else
                            ind_k(:,inds{k}) = true;
                        end
                        
                    elseif obj.dim==3
                        ind_k = false(grid_lengths(1),grid_lengths(2),grid_lengths(3));
                        
                        if obj.splitting_dim==1
                            ind_k(inds{k},:,:) = true;
                        elseif obj.splitting_dim==2
                            ind_k(:,inds{k},:) = true;
                        else
                            ind_k(:,:,inds{k}) = true;
                        end
                    end
                    
                    ind_k = ind_k(:);
                    
                    ind_k = in_index(ind_k);
                    
                    ii = [ii;ind_k(ii_k)];
                    jj = [jj;jj_k+step];
                    zz = [zz;zz_k];
                    
                end
                step = step + obj.children{k}.length();
            end
        end
        
        %  interpSparseMatrixZone_vecs(obj,X)
        %  Creates a sparse interpolating matrix along the zones
        %  for a list of points X.
        %
        %     Input:
        %      grid: cell array of grid values.
        %
        %    Output:
        %  ii,jj,kk: vectors used to construct the sparse matrix.
        function [ii,jj,zz,sum] = interpMatrixGrid_vecs(obj,grid)
            
            ii = []; jj = []; zz = [];
            grid_lengths = cellfun(@(x)length(x),grid);
            
            [n,~] = size(grid_lengths);
            
            sum = zeros(grid_lengths);
            
            if n>1
                grid_lengths = grid_lengths';
            end
            
            sub_grids = cell(2,1);
            
            inds = cell(2,1);
            
            not_empty_child = true(2,1);
            
            for k=1:2
                [sub_grids{k},inds{k}] = obj.children{k}.IndDomainGrid(grid);
                not_empty_child(k) = all(cellfun(@(x)~isempty(x),sub_grids{k}));
            end
            
            in_index = (1:prod(grid_lengths))';
            
            step = 0;
            
            for k=1:2
                if not_empty_child(k) && ~obj.children{k}.is_null
                    if obj.children{k}.is_leaf
                        
                        sumk = obj.children{k}.evalfGridBump(sub_grids{k});
                        [ii_k,jj_k,zz_k] = find(obj.children{k}.interpMatrixGrid(sub_grids{k}));
                        W = sumk(:);
                        
                        zz_k = zz_k.*W(ii_k);
                        
                        if size(ii_k,2)~=1
                            ii_k = ii_k';
                        end
                        
                        if size(jj_k,2)~=1
                            jj_k = jj_k';
                        end
                        
                        if size(zz_k,2)~=1
                            zz_k = zz_k';
                        end
                        
                    else
                        [ii_k,jj_k,zz_k,sumk] = obj.children{k}.interpMatrixGrid_vecs(sub_grids{k});
                    end
                    
                    if obj.dim==2
                        ind_k = false(grid_lengths(1),grid_lengths(2));
                        ind_k(inds{k}{1},inds{k}{2}) = true;
                        
                        sum(inds{k}{1},inds{k}{2}) = sum(inds{k}{1},inds{k}{2})+sumk;
                        
                    elseif obj.dim==3
                        ind_k = false(grid_lengths(1),grid_lengths(2),grid_lengths(3));
                        ind_k(inds{k}{1},inds{k}{2},inds{k}{3}) = true;
                        
                        sum(inds{k}{1},inds{k}{2},inds{k}{3}) = sum(inds{k}{1},inds{k}{2},inds{k}{3})+sumk;
                    end
                    
                    ind_k = ind_k(:);
                    
                    ind_k = in_index(ind_k);
                    
                    ii = [ii;ind_k(ii_k)];
                    jj = [jj;jj_k+step];
                    zz = [zz;zz_k];
                    
                end
                step = step + obj.children{k}.length();
            end
        end
        
        
        %  evalfZoneGrid(obj,X)
        %  Evaluates a grid along the zones for a grid X.
        %
        %     Input:
        %      grid: cell array of grid values.
        %
        %    Output:
        %      vals: values of approximation on grid along the zones.
        function vals = evalfZoneGrid(obj,X)
            grid_lengths = cellfun(@(x)length(x),X);
            
            [n,~] = size(grid_lengths);
            
            if n>1
                grid_lengths = grid_lengths';
            end
            
            %Put the order at the end; this makes things
            %easier to deal with in matlab
            vals = zeros(grid_lengths);
            
            sub_grids = cell(2,1);
            
            inds = cell(2,1);
            
            not_empty_child = true(2,1);
            
            sub_grids{1} = X;
            sub_grids{2} = X;
            
            midpoint = mean(obj.zone(obj.splitting_dim,:));
            
            inds{1} = X{obj.splitting_dim} < midpoint;
            inds{2} = ~inds{1};
            
            for k=1:2
                sub_grids{k}{obj.splitting_dim} = sub_grids{k}{obj.splitting_dim}(inds{k});
                not_empty_child(k) = all(cellfun(@(x)~isempty(x),sub_grids{k}));
            end
            
            
            %calculate values for the children
            for k=1:2
                if  not_empty_child(k) && ~obj.children{k}.is_null
                    if ~obj.children{k}.is_leaf
                        child_vals = obj.children{k}.evalfZoneGrid(sub_grids{k});
                    else
                        child_vals = obj.children{k}.evalfGrid(sub_grids{k});
                    end
                    
                    if obj.dim== 2
                        if obj.splitting_dim == 1
                            vals(inds{k},:) = child_vals;
                        else
                            vals(:,inds{k}) = child_vals;
                        end
                    elseif obj.dim == 3
                        if obj.splitting_dim == 1
                            vals(inds{k},:,:) = child_vals;
                        elseif obj.splitting_dim == 2
                            vals(:,inds{k},:) = child_vals;
                        else
                            vals(:,:,inds{k}) = child_vals;
                        end
                    end
                end
            end
        end
        
        %  interpMatrixZone_vecs(obj,X)
        %  Creates a sparse interpolating matrix for a list of points X.
        %
        %     Input:
        %      grid: cell array of grid values.
        %
        %    Output:
        %  ii,jj,kk: vectors used to construct the sparse matrix.
        function [ii,jj,zz] = interpMatrixZone_vecs(obj,X)
            
            [numpts,~] = size(X);
            
            ind = false(numpts,2);
            
            ii = [];
            jj = [];
            zz = [];
            
            mid_point = mean(obj.zone(obj.splitting_dim,:));
            
            ind(:,1) = X(:,obj.splitting_dim)<mid_point;
            ind(:,2) = ~ind(:,1);
            
            step = 0;
            
            %calculate values for the children
            for k=1:2
                if any(ind(:,k)) && ~obj.children{k}.is_null
                    
                    in_index = (1:numpts)';
                    in_index = in_index(ind(:,k));
                    
                    if obj.children{k}.is_leaf
                        M = obj.children{k}.interpMatrixPoints(X(ind(:,k),:));
                        [iik,jjk,zzk] = find(M);
                        
                        if(size(iik,2)~=1)
                            iik = iik.';
                            jjk = jjk.';
                            zzk = zzk.';
                        end
                        
                    else
                        [iik,jjk,zzk] = obj.children{k}.interpMatrixZone_vecs(X(ind(:,k),:));
                    end
                    
                    ii = [ii;in_index(iik)];
                    jj = [jj;jjk+step];
                    zz = [zz;zzk];
                end
                step = step+length(obj.children{k});
            end
        end
        
        %  interpMatrix_vecs(obj,X)
        %  Creates a sparse interpolating matrix for a list of points X.
        %
        %     Input:
        %      grid: cell array of grid values.
        %
        %    Output:
        %  ii,jj,kk: vectors used to construct the sparse matrix.
        function [ii,jj,zz,sum] = interpMatrix_vecs(obj,X)
            
            [numpts,~] = size(X);
            
            ind = false(numpts,2);
            
            sum = zeros(size(X,1),1);
            
            ii = [];
            jj = [];
            zz = [];
            
            ind(:,1) = obj.children{1}.InDomain(X);
            ind(:,2) = obj.children{2}.InDomain(X);
            
            step = 0;
            
            %calculate values for the children
            for k=1:2
                if any(ind(:,k)) && ~obj.children{k}.is_null
                    
                    in_index = (1:numpts)';
                    in_index = in_index(ind(:,k));
                    
                    if obj.children{k}.is_leaf
                        
                        sumk = obj.children{k}.evalfBump(X(ind(:,k),:));
                        M = obj.children{k}.interpMatrixPoints(X(ind(:,k),:));
                        [iik,jjk,zzk] = find(M);
                        zzk = zzk.*sumk(iik);
                        
                        if(size(iik,2)~=1)
                            iik = iik.';
                            jjk = jjk.';
                            zzk = zzk.';
                        end
                        
                    else
                        [iik,jjk,zzk,sumk] = obj.children{k}.interpMatrix_vecs(X(ind(:,k),:));
                    end
                    
                    ii = [ii;in_index(iik)];
                    jj = [jj;jjk+step];
                    zz = [zz;zzk];
                    sum(ind(:,k)) = sum(ind(:,k)) + sumk;
                end
                step = step+length(obj.children{k});
            end
        end
        
        %  evalfZone(obj,X)
        %  Evaluates the approximation along the zones for a list of
        %  points.
        %
        %     Input:
        %         X: list of points.
        %
        %    Output:
        %      vals: values of approximation along the zones
        %            for the list of points.
        function vals = evalfZone(obj,X)
            
            [numpts,~] = size(X);
            
            vals = zeros(numpts,1);
            
            ind = false(numpts,2);
            
            mid_point = mean(obj.zone(obj.splitting_dim,:));
            
            ind(:,1) = X(:,obj.splitting_dim)<mid_point;
            ind(:,2) = ~ind(:,1);
            
            child_vals = cell(2,1);
            
            %calculate values for the children
            for k=1:2
                if any(ind(:,k)) && ~obj.children{k}.is_null
                    if obj.children{k}.is_leaf
                        child_vals{k} = obj.children{k}.evalf(X(ind(:,k),:));
                    else
                        child_vals{k} = obj.children{k}.evalfZone(X(ind(:,k),:));
                    end
                else
                    child_vals{k} = [];
                end
            end
            
            vals(ind(:,1)) = child_vals{1};
            vals(ind(:,2)) = child_vals{2};
        end
        

        
%         function IsGeometricallyRefined = IsGeometricallyRefined(obj)
%             G1 = obj.children{1}.IsGeometricallyRefined();
%             G2 = obj.children{2}.IsGeometricallyRefined();
%             IsGeometricallyRefined = G1 & G2;
%         end
        
        % Plots the domains of the children.
        function plotdomain(obj)
            obj.children{1}.plotdomain();
            obj.children{2}.plotdomain();
        end
        
        % Plots the zones of the children.
        function plotzone(obj)
            if(~isempty(obj.children{1}.zone))
                obj.children{1}.plotzone();
            end
            
            if(~isempty(obj.children{2}.zone))
                obj.children{2}.plotzone();
            end
        end
        
        %Coarsens the leaves of the patch.
        function Coarsen(obj)
            for k=1:2
                obj.children{k}.Coarsen();
            end
            obj.cheb_length = length(obj.children{1})+length(obj.children{2});
            obj.iscoarse = obj.children{1}.iscoarse & obj.children{2}.iscoarse;
        end
        
        %Coarsens the leaves of the patch.
        function cvals = Fine2Coarse(obj,vals,j)
            for k=1:2
                cvals1 = obj.children{1}.Fine2Coarse(vals(1:length(obj.children{1})),j);
                cvals2 = obj.children{2}.Fine2Coarse(vals(length(obj.children{1})+1:end),j);
            end
            cvals = [cvals1;cvals2];
        end
        
        function rvals = Coarse2Fine(obj,vals)
            for k=1:2
                rvals1 = obj.children{1}.Coarse2Fine(vals(1:length(obj.children{1})));
                rvals2 = obj.children{2}.Coarse2Fine(vals(length(obj.children{1})+1:end));
            end
            rvals = [rvals1;rvals2];
        end
        
        %Refines the leaves of the patch.
        function Refine(obj)
            for k=1:2
                obj.children{k}.Refine();
            end
            obj.cheb_length = length(obj.children{1})+length(obj.children{2});
            obj.iscoarse = obj.children{1}.iscoarse & obj.children{2}.iscoarse;
        end
        
        %Method returns vector of values of interpolating points of patch.
        function vals = Getvalues(obj)
            vals = [obj.children{1}.Getvalues();obj.children{2}.Getvalues()];
        end
        
        %Recursive method that splits the children of a patch along a given
        %dimension.
        %
        %       Input:
        %   split_dim: splitting dimension
        %    set_vals: indicator if new children will have values
        %              interpolated from the parent.
        function split(obj,split_dim,set_vals)
            
            if nargin==2
                set_vals = false;
            end
            
            for k=1:2
                if obj.children{k}.is_leaf && ~obj.children{k}.is_null
                    obj.children{k} = obj.children{k}.split(split_dim,set_vals);
                elseif ~obj.children{k}.is_null
                    obj.children{k}.split(split_dim,set_vals);
                end
            end
            obj.cheb_length = length(obj.children{1})+length(obj.children{2});
            obj.domain = [obj.children{1}.domain(:,1) obj.children{2}.domain(:,2)];
        end
        
        % String function for method.
        function str = toString(obj)
            str = strvcat(strcat('1',obj.children{1}.toString()),strcat('2',obj.children{2}.toString()));
        end
        
        function reset(obj)
            obj.is_refined = false;
            obj.is_geometric_refined = false;
            
            if ~obj.is_leaf
                reset(obj.children{1});
                reset(obj.children{2});
            end
        end
        
        function Setvalues(obj,f)
            if isnumeric(f)
                obj.children{1}.Setvalues(f(1:length(obj.children{1})));
                obj.children{2}.Setvalues(f(length(obj.children{1})+1:end));
            else
                obj.children{1}.Setvalues(f);
                obj.children{2}.Setvalues(f);
            end
        end
        
    end
    methods(Access = protected)
        

        
        % Override copyElement method:
        function cpObj = copyElement(obj)
            % Make a shallow copy of all four properties
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            % Make a deep copy of the DeepCp object
            if ~obj.is_leaf
                cpObj.children{1} = obj.children{1}.copy();
                cpObj.children{2} = obj.children{2}.copy();
            end
        end
    end
    
end

