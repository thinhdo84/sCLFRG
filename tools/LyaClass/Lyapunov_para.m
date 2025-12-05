classdef Lyapunov_para< handle
    properties (SetAccess = private)
        V_list;
        alpha_list;
        N; %number of basis Lyapunov functions
    end
    methods (Access = public)
        function obj = Lyapunov_para(V_list,alpha_list)
            % Construct an instance of this class
            % Detailed explanation goes here
            obj.N = numel(V_list);
            obj.alpha_list = alpha_list;
            obj.V_list = V_list;
        end
        
        function val=eval(obj,X)
            val=0;
            for i=1:obj.N
                if obj.alpha_list(i)>1e-12
                    val = val + obj.alpha_list(i) * obj.V_list{i}.eval(X);
                end
            end
        end
        
        function dV=grad(obj,X)
            nx=numel(X);
            dV = zeros(1,nx);
            for i=1:obj.N
                if obj.alpha_list(i)>0
                    dV = dV + obj.alpha_list(i) * obj.V_list{i}.grad(X);
                end
            end
        end
    end
end

