%Class of tanh NN Lyapunov function
classdef Lyapunov_NNtanh< handle
    properties (SetAccess = private)
        Weight; bias; %weights and biases
        nx; nV; nHiddenLayer;
        L; % the size of the hidden layer
    end
    methods (Access = public)
        function obj = Lyapunov_NNtanh(Weight, bias)
            nx = size(Weight{1},2);
            
            nV = size(bias{end},1);
            nHiddenLayer = numel(Weight) - 1;
            if nHiddenLayer~=1
                error('this code is only for single hidden layer network')
            end
            % Construct an instance of this class
            % Detailed explanation goes here
            obj.nx = nx;
            obj.nV = nV;
            obj.nHiddenLayer = nHiddenLayer;
            obj.Weight = Weight;
            obj.bias = bias;
            obj.L = size(Weight{1},1);
        end
        
        function val=eval(obj,X)
            y1 = Vec_tanhActive(obj.Weight{1}*X +  obj.bias{1});
            val = obj.Weight{2}*y1 +  obj.bias{2};
        end

        
        function dV = grad(obj,X)
            dV = zeros(1,obj.nx);
            for k = 1:obj.nx
                for i =1:obj.L
                    dV(k) = dV(k)...
              +obj.Weight{2}(i)*dtanh(obj.Weight{1}(i,:)*X+obj.bias{1}(i))*obj.Weight{1}(i,k);
                end
            end
        end
    end
end


function y = tanhActive(x)
y = tanh(x);
end
function y = Vec_tanhActive(x)
if ~isvector(x)
    error('input must be a vector')
end
n = numel(x);
y = zeros(n,1);
for i = 1:n
    y(i) = tanhActive(x(i));
end
end

function dy = dtanh(x)
dy = 1-tanhActive(x)^2;
end
