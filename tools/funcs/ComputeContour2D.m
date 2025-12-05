function xb = ComputeContour2D(rBound,V,levelValue, varargin)
% Create input parser
p = inputParser;
% Add optional parameters with default values
addParameter(p, 'Nsample', 120);
% Parse inputs
parse(p, varargin{:});
% Get the values
Nsample = p.Results.Nsample;


tt = linspace(0,2*pi+0.05,Nsample);
X = rBound*[cos(tt); sin(tt)];
opts = optimoptions('fmincon', 'Display', 'none', ...
                    'Algorithm', 'interior-point', 'MaxFunctionEvaluations', 1e4); 
f_obj = @(lambda) -lambda;

xb = zeros(2,numel(tt));
for i = 1:numel(tt)
    constraint = @(v) deal([], V.eval(v*X(:,i))-levelValue);
    [lambda_star, ~] = fmincon(f_obj, 1e-4, [], [], [], [], 0, 1, constraint, opts);
    xb(:,i) = lambda_star*X(:,i);
end


end