% +linearSolvers/sor.m
function [x, iterations, final_error] = sor(A, b, omega, x0, tol, max_iter)
% Solves Ax = b using the Successive Over-Relaxation (SOR) method.

    n = size(A, 1);
    if nargin < 4 || isempty(x0)
        x = zeros(n, 1); % Initial guess
    else
        x = x0(:); % Ensure column vector
    end
    if nargin < 5
        tol = 1e-6;      % Default tolerance
    end
    if nargin < 6
        max_iter = 1000; % Default max iterations
    end

    % Check for zero diagonal elements
    if any(abs(diag(A)) < eps)
        error('Matrix has zero on the diagonal. SOR requires non-zero diagonal elements.');
    end
    
    D = diag(diag(A));
    L = tril(A, -1); % Strictly lower triangle
    U = triu(A, 1);  % Strictly upper triangle
    
    % Precompute inverse for matrix form (can be slow)
    % More efficient implementations solve the triangular system in each step
    inv_DLw = inv(D + omega * L);
    B_sor = inv_DLw * ((1 - omega) * D - omega * U);
    f_sor = omega * inv_DLw * b(:);

    x_prev = x;
    iterations = 0;
    final_error = inf;

    for k = 1:max_iter
        iterations = k;
        % Matrix form: x(k+1) = B_sor * x(k) + f_sor
        x = B_sor * x_prev + f_sor;

        % Check for convergence
        final_error = norm(x - x_prev, inf);
        if final_error < tol
            break;
        end
        x_prev = x;
    end

    if iterations == max_iter && final_error >= tol
        warning('SOR method (omega=%.2f) did not converge within %d iterations. Error = %e', omega, max_iter, final_error);
    end
end
