% function sigma = support_functional_smoothPolytope(H, p, y)
%     % Computes the support functional sigma_U(y)
%     % where U = {u | norm(H*u,2p)^(2p) <= 1}
%     y = y(:);
%     % Check dimensions
%     [r, m] = size(H);
%     assert(length(y) == m, 'y must be of compatible dimension with H');
%
%     % Check rank
%     if rank(H) < m
%         error('H must be full column rank');
%     end
%
%     % Compute q = Hölder conjugate of 2p
%     q = 2*p / (2*p - 1);
%
%     % Compute the linear map: H (H^T H)^{-1} y
%     H_pinv = H * ((H' * H)^(-1));  % equivalent to H * (H^T H)^{-1}
%     v = H_pinv * y;
%
%     % Compute q-norm of the result
%     sigma = norm(v, q);
% end


function sigma = support_functional_smoothPolytope(H, p, y)
% Computes the support function of the set U = {u : ||H u||_{2p}^{2p} <= 1}
% for square, invertible H
%
% Inputs:
%   H : [m x m] square matrix
%   y : [m x 1] direction vector
%   p : positive scalar, norm exponent (e.g., p = 5 for ||·||_{10})
%
% Output:
%   sigma : support function value in direction y
y = y(:);
[m, n] = size(H);

% Check that H is square
if m ~= n
    error('H must be square (m x m)');
end

% Check that H is invertible
if rank(H) < m
    error('H must be full rank (invertible)');
end

% Compute dual norm exponent
q = 2 * p / (2 * p - 1);
% Compute sigma
HTinv_y = (H^(-1))'*y        % (H^{-1})^T * y
sigma = norm(HTinv_y, q);  % Dual norm

end
