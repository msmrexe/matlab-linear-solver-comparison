% +linearSolvers/forwardSubstitution.m
function y = forwardSubstitution(L, b)
% Solves the lower triangular system Ly = b using forward substitution.

    n = size(L, 1);
    y = zeros(n, 1);
    b = b(:); % Ensure b is column vector

    if abs(L(1,1)) < eps
        error('Matrix L is singular (zero on diagonal).');
    end
    y(1) = b(1) / L(1,1);

    for i = 2:n
        if abs(L(i,i)) < eps
             error('Matrix L is singular (zero on diagonal).');
        end
        sum_ly = L(i, 1:i-1) * y(1:i-1); % Vectorized sum
        y(i) = (b(i) - sum_ly) / L(i,i);
    end
end
