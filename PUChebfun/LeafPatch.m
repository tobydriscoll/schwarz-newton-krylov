classdef LeafPatch<Patch
    % This is the abstract class for a leaf object. This is used with the
    % PUPatch object.
    
    % LeafPatch(var_struct) serves as the base constructor for a leafPatch
    % object, where var_struct is a structure with fields:
    %
    % *var_struct.domain, the domain used for the Chebyshev polynomial.
    %
    % *var_struct.zone, the zone (non overlapping part from partition). If the
    % zone is undefined, the zone is set to obj.domain.
    %
    % *var_struct.outerbox, the domain of the function. If the outerbox is
    % undefined, the outerbox is set to obj.domain.
    properties
        index = [];
        chebweights = [];
        %cheb_bump = chebfun( {0,@(x) exp(1-1./(1-x.^2)),0},[-20 -1 1 20]);
        bump
        values = [];
        coeffs
        is_interp
        % deg_in %index for the standard degrees
        degs %array of degrees along the dimensions
        cdegs %coarse degrees
        orig_degs %original max degrees
        orig_cdegs %coarse max degrees
        C2FinterpMat %Coarse to fine interpolation matrix
    end
    
    properties (Constant)
        %standard_variables = load('cheb_points_matrices.mat');
        %standard_degs = [3 5 9 17 33 65 129];
        %cheb_bump = @(x) exp(1-1./(1-x.^2));
        cheb_bump = chebfun( {0,@(x) exp(1-1./(1-x.^2)),0},[-20 -1 1 20]);
        invf = @(x,dom) 2/diff(dom)*x-sum(dom)/diff(dom); %takes points from a domain to [-1 1]
        forf = @(x,dom) 0.5*diff(dom)*x+0.5*sum(dom); %takes points from [-1 1] to a domain
    end
    
    methods
        function obj = LeafPatch(var_struct)
            
            if isfield(var_struct, 'domain')
                obj.domain = var_struct.domain;
                [obj.dim,~] = size(obj.domain);
            else
                error('A domain needs to be specified.');
            end
            
            obj.split_flag = true(obj.dim,1);
            obj.degs  = zeros(1,obj.dim);
            obj.cdegs  = zeros(1,obj.dim);
            obj.orig_degs  = zeros(1,obj.dim);
            obj.orig_cdegs  = zeros(1,obj.dim);
            
            obj.tol = 1e-12;
            obj.degs(:) = 128;
            obj.cdegs(:) = 9;
            
            obj.orig_degs(:) = 128;
            obj.orig_cdegs(:) = 9;
            
            if isfield(var_struct, 'outerbox')
                obj.outerbox = var_struct.outerbox;
            else
                obj.outerbox = obj.domain;
            end
            
            if isfield(var_struct, 'zone')
                obj.zone = var_struct.zone;
            else
                obj.zone = obj.domain;
            end
            
            if isfield(var_struct, 'degs')
                obj.degs = var_struct.degs;
            end
            
            if isfield(var_struct, 'cdegs')
                obj.cdegs = var_struct.cdegs;
            end
            
            if isfield(var_struct, 'orig_degs')
                obj.orig_degs = var_struct.orig_degs;
            end
            
            if isfield(var_struct, 'orig_cdegs')
                obj.orig_cdegs = var_struct.orig_cdegs;
            end
            
            if isfield(var_struct, 'split_flag')
                obj.split_flag = var_struct.split_flag;
            end
            
            if isfield(var_struct, 'tol')
                obj.tol = var_struct.tol;
            end
            

            
            
            
            obj.is_geometric_refined = true; %square is always geometrically refined
            [obj.dim,~] = size(obj.domain);
            obj.is_leaf = true;
            
            obj.bump = cell(3,1);
            
            %             for k=1:obj.dim
            %                 % if isequal(obj.domain(k,:),obj.outerbox(k,:))
            %                 %     w = @(x) ones(size(x));
            %                 if isequal(obj.domain(k,:),obj.outerbox(k,:))
            %                     w = @(x)obj.cheb_bump((obj.invf(x,obj.domain(k,:))+1)/2);
            %                     obj.bump{k} = @(x) (x<=obj.domain(k,1)) + (x>=obj.domain(k,2)) + (x>obj.domain(k,1) & x<obj.domain(k,2)).*w(x);
            %                 elseif obj.domain(k,1) == obj.outerbox(k,1)
            %                     w = @(x)obj.cheb_bump((obj.invf(x,obj.domain(k,:))+1)/2);
            %                     obj.bump{k} = @(x) (x<=obj.domain(k,1)) + (x>obj.domain(k,1) & x<obj.domain(k,2)).*w(x);
            %                 elseif obj.domain(k,2) == obj.outerbox(k,2)
            %                     w = @(x) obj.cheb_bump((obj.invf(x,obj.domain(k,:))-1)/2);
            %                     obj.bump{k} = @(x)(x>=obj.domain(k,2))+(x>obj.domain(k,1) & x<obj.domain(k,2)).*w(x);
            %                 else
            %                     w = @(x) obj.cheb_bump(obj.invf(x,obj.domain(k,:)));
            %                     obj.bump{k} = @(x)(x>obj.domain(k,1) & x<obj.domain(k,2)).*w(x);
            %                 end
            %             end
            
            
            for k=1:obj.dim
                % if isequal(obj.domain(k,:),obj.outerbox(k,:))
                %     w = @(x) ones(size(x));
                if obj.domain(k,1) == obj.outerbox(k,1)
                    w = @(x) obj.cheb_bump((obj.invf(x,obj.domain(k,:))+1)/2);
                elseif obj.domain(k,2) == obj.outerbox(k,2)
                    w = @(x) obj.cheb_bump((obj.invf(x,obj.domain(k,:))-1)/2);
                else
                    w = @(x) obj.cheb_bump(obj.invf(x,obj.domain(k,:)));
                end
                obj.bump{k} = w;
            end
        end
        
        
        % Evaluates the bump function of a patch.
        %
        %  Input:
        %      X: set of points to evaluate at
        %
        % Output:
        %     ef: array of length(X) of the patch's weight evaluated at X.
        function ef = evalfBump(obj,X)
            
            ef = ones(size(X,1),1);
            for k=1:obj.dim
                ef = ef.*obj.bump{k}(X(:,k));
            end
            
        end
        
        % Evaluates the bump function of a patch.
        %
        %  Input:
        %      X: cellarray of grids to be evaluated on.
        %
        % Output:
        %     ef: array of dim(X) of the patch's weight evaluated at X.
        function ef = evalfGridBump(obj,X)
            
            W = cell(3,1);
            
            for i=1:obj.dim
                W{i} = obj.bump{i}(X{i});
            end
            
            if obj.dim==1
                ef = W{1};
            elseif obj.dim==2
                ef = W{1}*W{2}.';
            else
                ef = reshape(W{3},1,1,length(W{3})).*(W{2}'.*W{1});
            end
            
        end
        
        % Plots the overlaping domain of the tree. Works only
        % in 2D and 3D.
        %
        %     Input:
        %     color: color of marker for the center of the domain.
        %
        function plotdomain(obj,color)
            
            if nargin==1
                color = 'black';
            end
            
            if obj.dim==2
                hold on;
                lengths = [diff(obj.domain(1,:));diff(obj.domain(2,:))];
                rectangle('position',[obj.domain(:,1)' lengths'],'LineWidth',2,'EdgeColor',color);
                plot(mean(obj.domain(1,:)),mean(obj.domain(2,:)),'.','MarkerSize',10,'Color',color);
                hold off;
            elseif obj.dim==3
                hold on;
                lengths = [diff(obj.domain(1,:));diff(obj.domain(2,:));diff(obj.domain(3,:))];
                center = sum(obj.domain,2)/2;
                %Vertices for Line Cube. Order matters
                X = [-1 -1 1 1 -1 -1 1 1 1 1 1 1 -1 -1 -1 -1 -1]';
                Y = [-1 1 1 -1 -1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1]';
                Z = [-1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 1 -1 1 1 -1]';
                %Example two cube matrix. Unit cube and one scaled/translated cube
                X1 = X*lengths(1)/2+center(1);
                Y1 = Y*lengths(2)/2+center(2);
                Z1 = Z*lengths(3)/2+center(3);
                %Single plot command for all 'cube lines'
                plot3(X1,Y1,Z1,'color','black');
                plot3(mean(obj.domain(1,:)),mean(obj.domain(2,:)),mean(obj.domain(3,:)),'.','MarkerSize',10,'Color','black');
                hold off;
            end
        end
        
        
        % Plots the zones of the tree. Works only
        % in 2D and 3D.
        %
        function plotzone(obj)
            
            if obj.dim==2
                hold on;
                lengths = [diff(obj.zone(1,:));diff(obj.zone(2,:))];
                rectangle('position',[obj.zone(:,1)' lengths'],'LineWidth',2);
                hold off;
            elseif obj.dim==3
                hold on;
                lengths = [diff(obj.zone(1,:));diff(obj.zone(2,:));diff(obj.zone(3,:))];
                center = sum(obj.zone,2)/2;
                %Vertices for Line Cube. Order matters
                X = [-1 -1 1 1 -1 -1 1 1 1 1 1 1 -1 -1 -1 -1 -1]';
                Y = [-1 1 1 -1 -1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1]';
                Z = [-1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 1 -1 1 1 -1]';
                %Example two cube matrix. Unit cube and one scaled/translated cube
                X1 = X*lengths(1)/2+center(1);
                Y1 = Y*lengths(2)/2+center(2);
                Z1 = Z*lengths(3)/2+center(3);
                %Single plot command for all 'cube lines'
                plot3(X1,Y1,Z1,'color','black');
                hold off;
            end
        end
        
        % Returns the points of the function
        function pts = points(obj)
            
            C = cell(obj.dim,1);
            
            for i=1:obj.dim
                C{i} = chebpts(obj.degs(i),obj.domain(i,:));
            end
            [out{1:obj.dim}] = ndgrid(C{:});
            
            pts = zeros(numel(out{1}),obj.dim);
            
            for i=1:obj.dim
                pts(:,i) = out{i}(:);
            end
            
        end
        
        % Returns a cell array of the grids of the domain
        function grid = leafGrids(obj)
            grid = cell(1,obj.dim);
            
            for i=1:obj.dim
                grid{i} = chebpts(obj.degs(i),obj.domain(i,:));
            end
        end
        
        function ln=length(obj)
            ln = prod(obj.degs);
        end
        
        % TODO. Figure out what to do here!
        function ef = evalf(obj,X,G)
            if nargin<3
                G = obj.coeffs;
            end
            
            [num_pts,~] = size(X);
            
            ef = zeros(num_pts,1);
            
            for i=1:num_pts
                ef(i) = evalfGrid(obj,num2cell(X(i,:)),G);
            end
            
        end
        
        function ef = evalfGrid(obj,X,G)
            
            if nargin<3
                G = obj.coeffs;
            end
            
            [~,d_ind] = sort(obj.degs,'descend');
            
            for k=1:obj.dim
                %Shift the points to the domain [-1 1]x[-1 1]
                X{d_ind(k)} = obj.invf(X{d_ind(k)},obj.domain(d_ind(k),:));
                
                %Evaluate the points at the Chebyshev polynomials
                F = chebtech.clenshaw(X{d_ind(k)},eye(obj.degs(d_ind(k))));
                
                %Multiply the coefficients with F
                G = chebfun3t.txm(G, F, d_ind(k));
            end
            
            ef  = G;
            
        end
        
        function C = DiffCoeffs(obj,diff_dim,order)
            
            if nargin<3
                order = 1;
            end
            
            unContractedModes = [1:diff_dim-1, diff_dim+1:obj.dim];
            
            C = chebfun3t.unfold(obj.coeffs,diff_dim);
            
            for i=1:order
                C = obj.computeDiffCoeffs(C);
            end
            
            C = chebfun3t.fold(C,obj.degs,diff_dim,unContractedModes);
            
        end
        
        % Evaluates the derivative on a grid.
        %
        %  Input:
        %      diff_dim: dimension derivative is taken in
        %         order: order of derivative (up to 2 right now)
        %             X: cellarray of grids to be evaluated on.
        %
        % Output:
        %     ef: matrix of dim(X) containing the interpolated derivative values
        function ef = evalfDiffGrid(obj,diff_dim,order,X)
            if nargin<3
                order = 1;
            end
            
            unContractedModes = [1:diff_dim-1, diff_dim+1:obj.dim];
            
            G = chebfun3t.unfold(obj.coeffs,diff_dim);
            
            for i=1:order
                G = obj.computeDiffCoeffs(G);
            end
            
            G = chebfun3t.fold(G,obj.degs,diff_dim,unContractedModes)/diff(obj.domain(diff_dim,:)/2)^order;
            
            if nargin<4
                ef = G;
            else
                ef = evalfGrid(obj,X,G);
            end
        end
        
        % Evaluates the approximant and its derivatives.
        %
        %  Input:
        %      X: set of points to evaluate at
        %
        % Output:
        %     ef: length(X) array containing the interpolated
        function ef = Diff(obj,diff_dim,order,X)
            if nargin<3
                order = 1;
            end
            
            unContractedModes = [1:diff_dim-1, diff_dim+1:obj.dim];
            G = chebfun3t.unfold(obj.coeffs,diff_dim);
            
            for i=1:order
                G = obj.computeDiffCoeffs(G);
            end
            
            G = chebfun3t.fold(G,obj.degs,diff_dim,unContractedModes)/diff(obj.domain(diff_dim,:)/2)^order;
            
            if nargin<4
                ef = G;
            else
                ef = evalfDiffGrid(obj,diff_dim,order,X,G);
            end
        end
        
        %  interpMatrixPoints(obj,X)
        %  This method creates a interpolating matrix given a list of
        %  points.
        %
        %  Input:
        %      X: list of points.
        %
        % Output:
        %      M: interpolating matrix.
        function M = interpMatrixPoints(obj,X)
            
            M = zeros(size(X,1),length(obj));
            
            for i=1:size(X,1)
                M(i,:) = obj.interpMatrixGrid(num2cell(X(i,:)));
            end
            
        end
        
        %  interpMatrixGrid(obj,grid)
        %  This method creates a interpolating matrix given a grid.
        %
        %  Input:
        %      X: cell array of grid values.
        %
        % Output:
        %      M: interpolating matrix.
        function M = interpMatrixGrid(obj,grid)
            G = obj.leafGrids();
            if obj.dim==1
                if iscell(grid)
                    M = barymat(grid{1},G{1});
                else
                    M = barymat(grid,G{1});
                end
            elseif obj.dim==2
                M = kron(barymat(grid{2},G{2}),barymat(grid{1},G{1}));
            elseif obj.dim==3
                M = kron(barymat(grid{3},G{3}),kron(barymat(grid{2},G{2}),barymat(grid{1},G{1})));
            end
        end
        
        % Evaluates the approximant on a grid.
        %
        %  Input:
        %      X: cellarray of grids to be evaluated on.
        %
        % Output:
        %     ef: matrix of dim(X) containing the interpolated values
        function ef = evalfBaryGrid(obj,X,G)
            
            for k=1:obj.dim
                %Shift the points to the right domain
                points = obj.invf(X{k},obj.domain(k,:));
                %x = chebpts(obj.degs(k));
               % w = chebtech2.barywts(obj.degs(k));
                f = bary(points,eye(obj.degs(k)));
                
                if length(X{k})==1
                    f = f.';
                end
                
                G = chebfun3t.txm(G, f, k);
            end
            ef = G;
        end
        
        % Evaluates the approximant and its derivatives.
        %
        %  Input:
        %      X: set of points to evaluate at
        %
        % Output:
        %     ef: length(X) array containing the interpolated
        function ef = evalfBarry(obj,X,G)
            
            [num_pts,~] = size(X);
            
            ef = zeros(num_pts,1);
            
            for i=1:num_pts
                ef(i) = evalfGrid(obj,num2cell(X(i,:)),G);
            end
            
        end
        
        function Coarsen(obj)
            
            if ~obj.iscoarse
                
                obj.iscoarse = true;
                
                obj.swap_degs = obj.degs;
                
%                 obj.degs = obj.cdegs;
%                 
%                 grid = obj.leafGrids();
%                 
%                 obj.degs = obj.swap_degs;
%                 
%                 obj.values = obj.evalfBaryGrid(grid,obj.values);
                
                obj.degs = obj.cdegs;
                
                obj.cheb_length = prod(obj.cdegs);
            end
            
        end
        
        function cvals = Fine2Coarse(obj,vals,k)
            
            if 2==nargin
                k=0;
            end
            
            if ~obj.iscoarse
                
                if k==0
                    obj.swap_degs = obj.degs;
                    
                    obj.degs = obj.cdegs;
                    
                    grid = obj.leafGrids();
                    
                    obj.degs = obj.swap_degs;
                    
                    cvals = obj.evalfBaryGrid(grid,reshape(vals,obj.degs));
                    
                else
                    
                    vals = reshape(vals,obj.degs);
                    
                    if obj.dim==1
                        cvals = vals(1:2^k:end);
                    elseif obj.dim==2
                        cvals = vals(1:2^k:end,1:2^k:end);
                    else
                        cvals = vals(1:2^k:end,1:2^k:end,1:2^k:end);
                    end
                    
                end
                
                cvals = cvals(:);
            end
            
        end
        
        
        % Construct for the ChebPatch
        %
        % This method refines the patch, resampling the patches values
        % on a fine grid.
        function Refine(obj)
            if obj.iscoarse
                
                obj.iscoarse = false;
                
%                 obj.degs = obj.swap_degs;
%                 
%                 grid = obj.leafGrids();
%                 
%                 obj.degs = obj.cdegs;
%                 
%                 obj.values = obj.evalfBaryGrid(grid,obj.values);
                
                obj.degs = obj.swap_degs;
                
                obj.cheb_length = prod(obj.degs);
            end
        end
        
        function rvals = Coarse2Fine(obj,vals)
            
            if obj.iscoarse
%                 rvals = obj.C2FinterpMat*vals;
                 
                vals = reshape(vals,obj.degs);
                scoeffs = zeros(obj.swap_degs);

                U = chebtech2.vals2coeffs( vals );  
                U = chebtech2.vals2coeffs( U.' ).'; 
                
                scoeffs(1:obj.degs(1),1:obj.degs(2)) = U;
                
                
                rvals = chebtech2.coeffs2vals( scoeffs );
                rvals = chebtech2.coeffs2vals( rvals.' ).'; 
                
                rvals = rvals(:);
%                 

                
%                 obj.degs = obj.swap_degs;
%                 
%                 grid = obj.leafGrids();
%                 
%                 obj.degs = obj.cdegs;
%                 
%                 rvals = obj.evalfBaryGrid(grid,reshape(vals,obj.degs));
%                 
%                 rvals = rvals(:);
            end
            
        end
        
        
        function Setvalues(obj,f)
            if isnumeric(f)
                obj.values = reshape(f,obj.degs);
            else
                if obj.dim==1
                    obj.values = f(obj.points());
                elseif obj.dim==2
                    points = obj.points();
                    obj.values = reshape(f(points(:,1),points(:,2)),obj.degs);
                else
                    points = obj.points();
                    obj.values = reshape(f(points(:,1),points(:,2),points(:,3)),obj.degs);
                end
            end
        end
        
    end
    
    methods (Abstract)
        %This method will split the child, creating a new PUPatch. If the
        %obj does not need to split, the method returns obj.
        Child = splitleaf(obj);
    end
    
    methods (Static)
        function cout = computeDiffCoeffs(c)
            %COMPUTEDERCOEFFS   Recurrence relation for coefficients of derivative.
            %   C is the matrix of Chebyshev coefficients of a (possibly array-valued)
            %   CHEBTECH object.  COUT is the matrix of coefficients for a CHEBTECH object
            %   whose columns are the derivatives of those of the original.
            
            [n, m] = size(c);
            cout = zeros(n, m);                        % Initialize vector {c_r}
            w = repmat(2*(1:n-1)', 1, m);
            v = w.*c(2:end,:);                           % Temporal vector
            cout(n-1:-2:1,:) = cumsum(v(n-1:-2:1,:), 1); % Compute c_{n-2}, c_{n-4}, ...
            cout(n-2:-2:1,:) = cumsum(v(n-2:-2:1,:), 1); % Compute c_{n-3}, c_{n-5}, ...
            cout(1,:) = .5*cout(1,:);                    % Adjust the value for c_0
        end
    end
end
