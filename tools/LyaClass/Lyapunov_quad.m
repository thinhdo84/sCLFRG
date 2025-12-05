classdef Lyapunov_quad < handle
    % Lyapunov_quad Class of quadratic Lyapunov function
    properties (SetAccess = private)
        P;nx;nc;
    end   
    methods (Access = public)
        function obj = Lyapunov_quad(P)
            nx = size(P,1);
            % Construct an instance of this class
            % Detailed explanation goes here
            obj.nx = nx;
            obj.P = P;
        end
        
        function val=eval(obj,X)
            val=X'*obj.P*X;
        end
        
        function dV = grad(obj,X) 
%             X
           dV = 2*(obj.P*X)';%?
           
        end
    end
end

