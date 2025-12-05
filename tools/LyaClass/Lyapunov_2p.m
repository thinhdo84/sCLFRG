%Class of 2p-norm Lyapunov function
classdef Lyapunov_2p< handle
    properties (SetAccess = private)
        F;p;nx;nc;
    end   
    methods (Access = public)
        function obj = Lyapunov_2p(F,p)
            [nc,nx] = size(F);
            % Construct an instance of this class
            % Detailed explanation goes here
            obj.nx = nx;
            obj.nc = nc;
            obj.F = F;
            obj.p = p;
        end
        
        function val=eval(obj,X)
            val=0;
            for i=1:obj.nc
               val = val + (obj.F(i,:)*X)^(2*obj.p);
            end
        end
        
        function dV = grad(obj,X) 
           dV = zeros(1,obj.nx);
           for k =1:obj.nx
               for i=1:obj.nc
                    dV(k) = dV(k)+ 2*obj.p*obj.F(i,k)*(obj.F(i,:)*X)^(2*obj.p-1);
               end               
           end
        end
    end
end

