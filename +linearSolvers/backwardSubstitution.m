% +linearSolvers/backwardSubstitution.m
function x = backwardSubstitution(U, y)
% Solves the upper triangular system Ux = y using backward substitution.

    n = size(U, 1);
    x = zeros(n, 1);
    y = y(:); % Ensure y is column vector

    if abs(U(n,n)) < eps
         error('Matrix U is singular (zero on diagonal).');
    end
    x(n) = y(n) / U(n,n);

    for i = n-1:-1:1
        if abs(U(i,i)) < eps
             error('Matrix U is singular (zero on diagonal).');
        end
        sum_ux = U(i, i+1:n) * x(i+1:n); % Vectorized sum
        x(i) = (y(i) - sum_ux) / U(i,i);
    end
end
