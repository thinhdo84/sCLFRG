function xb = SamplePolytopeBoundary(rBound,A,b,varargin)
% Sample Nsample points from the polytope of {x: Ax<=b}
n = size(A,2);

% Create input parser
p = inputParser;
% Add optional parameters with default values
addParameter(p, 'Nsample', 120);
% Parse inputs
parse(p, varargin{:});
% Get the values
Nsample = p.Results.Nsample;


if n == 2
    
    tt = linspace(0,2*pi+0.05,Nsample);
    X = rBound*[cos(tt); sin(tt)];
    opts = optimoptions('fmincon', 'Display', 'none', ...
        'Algorithm', 'interior-point', 'MaxFunctionEvaluations', 1e4);
    f_obj = @(lambda) -lambda;
    
    xb = zeros(2,numel(tt));
    for i = 1:numel(tt)
        [lambda_star, ~] = fmincon(f_obj, 1e-4, A*X(:,i), b, [], [], 0, 1, [], opts);
        xb(:,i) = lambda_star*X(:,i);
    end
    
else
    error('n>2 is under dev')
end
end