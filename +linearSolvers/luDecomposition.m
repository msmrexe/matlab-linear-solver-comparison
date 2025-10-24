% +linearSolvers/luDecomposition.m
function x = luDecomposition(A, b)
% Solves Ax = b using LU Decomposition with partial pivoting (LUP).

    [n, ~] = size(A);
    b = b(:); % Ensure b is a column vector
    
    % Perform LU decomposition with pivoting using MATLAB's built-in function
    % [L, U, P] = lu(A) gives P*A = L*U
    [L, U, P] = lu(A);

    % Check if U is singular (matrix A is singular)
    if any(abs(diag(U)) < eps)
        % Check for consistency: P*A*x = P*b => L*U*x = P*b
        % If rank(L*U) < rank([L*U, P*b]), no solution
        % If rank(L*U) == rank([L*U, P*b]) < n, infinite solutions
        % Simplified check: Check if P*b lies in the column space of L*U
        % However, just erroring out for singular is common.
        error('Matrix is singular or nearly singular. Cannot find unique solution using LU.');
    end

    % Solve Ly = Pb using forward substitution
    y = L \ (P * b);
    
    % Solve Ux = y using backward substitution
    x = U \ y;

end
