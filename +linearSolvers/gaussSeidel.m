% +linearSolvers/gaussSeidel.m
function [x, iterations, final_error] = gaussSeidel(A, b, x0, tol, max_iter)
% Solves Ax = b using the Gauss-Seidel iterative method (matrix form).

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

    % Check for zero diagonal elements
    if any(abs(diag(A)) < eps)
        error('Matrix has zero on the diagonal. Gauss-Seidel requires non-zero diagonal elements for this implementation.');
    end

    D = diag(diag(A));
    L = tril(A, -1); % Strictly lower triangle
    U = triu(A, 1);  % Strictly upper triangle
    
    DL_inv = inv(D + L); % Precompute inverse (can be slow for large matrices)
                         % More efficient: use forward substitution in loop
                         % But for clarity/simplicity here, precompute.
    
    x_prev = x;
    iterations = 0;
    final_error = inf;

    for k = 1:max_iter
        iterations = k;
        % Matrix form: x(k+1) = (D+L)^-1 * (b - U*x(k))
        x = DL_inv * (b(:) - U * x_prev);
        
        % Check for convergence
        final_error = norm(x - x_prev, inf);
        if final_error < tol
            break;
        end
        x_prev = x;
    end
    
    if iterations == max_iter && final_error >= tol
       warning('Gauss-Seidel method did not converge within %d iterations. Error = %e', max_iter, final_error);
    end
end
