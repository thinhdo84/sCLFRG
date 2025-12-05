function [a,b] = tangent_line_params(c, x)
% c, x: 2x1 (or 1x2) vectors in R^2
% Returns a (2x1) and scalar b for the line a'*y = b
% that passes through x and is tangent to the circle centered at c
% with radius ||x-c||, and satisfies a'*c <= b.

c = c(:); x = x(:);
a = x - c;
if all(a == 0)
    error('x must differ from c (nonzero radius).');
end
b = a.' * x;
end
