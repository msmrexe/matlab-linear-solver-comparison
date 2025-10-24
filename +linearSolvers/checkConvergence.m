% +linearSolvers/checkConvergence.m
function can_apply = checkConvergence(A, method, omega)
% Checks if an iterative method is likely to converge for matrix A.
% Calculates the spectral radius of the iteration matrix.

    can_apply = false; % Default
    try
        n = size(A, 1);
        D = diag(diag(A));
        
        % Check for zero diagonal elements, which cause failure
        if any(abs(diag(D)) < eps)
             warning('Matrix has near-zero diagonal elements. Iterative methods likely unstable/fail.');
             return;
        end

        L = tril(A, -1);
        U = triu(A, 1);

        iteration_matrix = [];

        if strcmpi(method, 'Jacobi')
            D_inv = diag(1./diag(D));
            iteration_matrix = -D_inv * (L + U);
        elseif strcmpi(method, 'GaussSeidel')
            DL_inv = inv(D + L); % Warning: inverse can be slow/unstable
            iteration_matrix = -DL_inv * U;
        elseif strcmpi(method, 'SOR')
            if nargin < 3
                error('Omega value required for SOR convergence check.');
            end
            inv_DLw = inv(D + omega * L); % Warning: inverse can be slow/unstable
            iteration_matrix = inv_DLw * ((1 - omega) * D - omega * U);
        else
            error('Unknown method specified for convergence check.');
        end

        % Calculate spectral radius (maximum absolute eigenvalue)
        eigenvalues = eig(iteration_matrix);
        spectral_radius = max(abs(eigenvalues));

        if spectral_radius < 1
            can_apply = true;
        end

    catch ME
        warning('Convergence check failed: %s', ME.message);
        % Could fail due to singularity, etc.
        can_apply = false;
    end
end
