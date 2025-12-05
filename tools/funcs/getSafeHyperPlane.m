function [a_safe, b_safe] = getSafeHyperPlane(x,PObstacle, cObstacle) 
if numel(x) ~=2
    error('the code works only in 2D')
end
[~, xtouch, ~]= find_largest_touching_level_set_newton(cObstacle, PObstacle, x, eye(2));
[a_safe,b_safe] = tangent_line_params(x, xtouch);
end