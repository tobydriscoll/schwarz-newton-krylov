classdef PUWeights
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        weights
        diffweights
        diff2weights
        chebweights
    end
    
    properties (Constant)
        overlap = 0.1;
    end
    
    methods
        function obj = PUWeights()
            
            t = obj.overlap;
            
            bump = chebfun( {0,@(x) exp(1-1./(1-x.^2)),0},[-20 -1 1 20]);
            
            R = 1 + obj.overlap;
            bl = @(x) bump((x+1)/R);
            br = @(x) bump((x-1)/R);
            
            %These could be precalculated. This is what I will do.
           weight1 = chebfun({@(x) bl(x)./(bl(x)+br(x)),0},[-1 obj.overlap 1],'vectorize');
            
           % weight1 = chebfun({1,@(x)-0.5/t*x+0.5,0},[-1 -t t 1]);
            
            obj.chebweights = [weight1 1-weight1];
            
            obj.weights = {@(x) 1./(1+exp(4*(1+t)^2*x./((t-x).*(2+t-x).*(t+x).*(2+t+x)))), ...
                           @(x) 1./(1+exp(-4*(1+t)^2*x./((t-x).*(2+t-x).*(t+x).*(2+t+x))))};
            

            %obj.diffweights = @(x) (1+t).^2.*(t+(-1).*x).^(-2).*(2+t+(-1).*x).^(-2).*(t+x).^(-2).*(2+t+x).^(-2).*...
            %                        (t.^2.*(2+t).^2+2.*(2+t.*(2+t)).*x.^2+(-3).*x.^4).*...
            %                        sech(2.*(1+t).^2.*(t+(-1).*x).^(-1).*(2+t+(-1).*x).^(-1).*x.*(t+x).^(-1).*(2+t+x).^(-1)).^2;
            
            obj.diffweights = diff(obj.chebweights);
           obj.diff2weights = diff(obj.chebweights,2);
            
            %obj.diff2weights = @(x) (1+t).^2.*(t+(-1).*x).^(-3).*(2+t+(-1).*x).^(-3).*(t+x).^(-3).*(2+t+x).^(-3).*...
            %                        sech(2.*(1+t).^2.*(t+(-1).*x).^(-1).*(2+t+(-1).*x).^(-1).*x.*(t+x).^(-1).*(2+t+x).^(-1)).^2.*...
            %                        (4.*x.*(3.*t.^2.*(2+t).^2.*(2+t.*(2+t))+(8+(-1).*t.*(2+t).*((-8)+3.*t.*(2+t))).*x.^2+(-3).*(2+t.*(2+t)).*x.^4+3.*x.^6)...
            %                        +(-4).*(1+t).^2.*(t+(-1).*x).^(-1).*(2+t+(-1).*x).^(-1).*(t+x).^(-1).*(2+t+x).^(-1).*(t.^2.*(2+t).^2+2.*(2+t.*(2+t)).*x.^2+(-3).*x.^4).^2.*...
            %                        tanh(2.*(1+t).^2.*(t+(-1).*x).^(-1).*(2+t+(-1).*x).^(-1).*x.*(t+x).^(-1).*(2+t+x).^(-1)));
        end
        
        function ef = evalf(obj,X,MIDINV,k,split_dim,diff_j)
            
            t = obj.overlap;
            
            if nargin == 4
                split_dim =1;
                diff_j =0;
            end
            
            if nargin == 5
                diff_j=0;
            end
            
           %h = @(x) 2/(MIDINV(2)-MIDINV(1))*x-(MIDINV(2)+MIDINV(1))/(MIDINV(2)-MIDINV(1));
           %SCALE = 2/(MIDINV(2)-MIDINV(1));
           % h = @(x) 2*t/(MIDINV(2)-MIDINV(1))*x-t*(MIDINV(2)+MIDINV(1))/(MIDINV(2)-MIDINV(1));
           % SCALE = 2*t/(MIDINV(2)-MIDINV(1));
           h = @(x) 2/diff(MIDINV)*obj.overlap*x-sum(MIDINV)/diff(MIDINV)*obj.overlap;
           SCALE = sum(MIDINV)/diff(MIDINV)*obj.overlap;
           
            %collect points along the splitting dimension.
            x = X(:,split_dim);
            
            [num_pts,~] = size(X);
            
            ef = zeros(num_pts,1);
            
            X_CENTER = h(x(x>MIDINV(1) & x<MIDINV(2)));
            
            if diff_j==0
                    ef(x<=MIDINV(1)) = k==1;
                    
                    ef(x>MIDINV(1) & x<MIDINV(2)) = feval(obj.chebweights(:,k),X_CENTER);
                    
                    ef(x>=MIDINV(2)) = k==2;
                    
            elseif diff_j==1
                    ef(x>MIDINV(1) & x<MIDINV(2)) = SCALE*feval(obj.diffweights(:,k),X_CENTER);
            elseif diff_j==2
                    ef(x>MIDINV(1) & x<MIDINV(2)) = SCALE^2*feval(obj.diff2weights(:,k),X_CENTER);
            else
                    ef(x>MIDINV(1) & x<MIDINV(2)) = SCALE^diff_j*feval(diff(obj.weights(:,k),diff_j),X_CENTER);
            end
        end
    end
end

