% +linearSolvers/jacobi.m
function [x, iterations, final_error] = jacobi(A, b, x0, tol, max_iter)
% Solves Ax = b using the Jacobi iterative method (matrix form).

    n = size(A, 1);
    if nargin < 3 || isempty(x0)
        x = zeros(n, 1); % Initial guess
    else
        x = x0(:); % Ensure column vector
    end
    if nargin < 4
        tol = 1e-6;      % Default tolerance
    end
    if nargin < 5
        max_iter = 1000; % Default max iterations
    end

    D = diag(diag(A));
    LU = A - D; % L + U combined

    % Check for zero diagonal elements
    if any(abs(diag(D)) < eps)
        error('Matrix has zero on the diagonal. Jacobi requires non-zero diagonal elements.');
    end

    D_inv = diag(1./diag(D));
    
    x_prev = x;
    iterations = 0;
    final_error = inf;

    for k = 1:max_iter
        iterations = k;
        x = D_inv * (b(:) - LU * x_prev);

        % Check for convergence
        final_error = norm(x - x_prev, inf); % Use infinity norm
        if final_error < tol
            break;
        end
        x_prev = x;
    end

    if iterations == max_iter && final_error >= tol
        warning('Jacobi method did not converge within %d iterations. Error = %e', max_iter, final_error);
    end
end
