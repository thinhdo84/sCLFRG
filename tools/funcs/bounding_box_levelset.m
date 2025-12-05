function [ub, lb] = bounding_box_levelset(V, c, alpha)
% Computes axis-aligned bounding box of the level set:
%   {x : V(x) <= alpha}
% using fmincon.

% Output: ub, vb describing the set:
% {x: lb <= x <= ub};

n = V.V_list{end}.nx;
% Options for fmincon solver
opts = optimoptions('fmincon', 'Display', 'none', ...
                    'Algorithm', 'interior-point', 'MaxFunctionEvaluations', 1e4); 
% Initialization
ub = zeros(n,1); % upper bounds
lb = zeros(n,1); % lower bounds

% Solve for the furthest hyperplane (linear cost, 1 convex constraint)
for i = 1:n
    % Maximize v_i: equivalent to minimizing -v_i
    f_obj = @(v) -v(i);
    constraint = @(v) deal([], V.eval(v) - alpha);
    [v_star, ~] = fmincon(f_obj, c, [], [], [], [], [], [], constraint, opts);
    ub(i) = c(i) + v_star(i);

    % Minimizing v_i
    f_obj = @(v) v(i);
    [v_star, ~] = fmincon(f_obj, c, [], [], [], [], [], [], constraint, opts);
    lb(i) = c(i) + v_star(i);
end
end