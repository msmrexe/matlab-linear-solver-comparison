% +linearSolvers/luDecomposition.m
function x = luDecomposition(A, b)
% Solves Ax = b using from-scratch LU Decomposition with partial
% pivoting (LUP), followed by from-scratch forward and backward substitution.

    [n, ~] = size(A);
    b = b(:); % Ensure b is a column vector
    
    % Perform LUP factorization using the custom function
    % PA = LU
    try
        [L, U, P] = linearSolvers.lupFactorization(A);
    catch ME
        % Propagate error from factorization (e.g., singularity)
        rethrow(ME);
    end

    % Solve Ly = Pb using forward substitution
    Pb = P * b;
    try
        y = linearSolvers.forwardSubstitution(L, Pb);
    catch ME
        % Check for inconsistency potentially revealed here
         if contains(ME.message, 'singular')
             % If L is singular but P*A wasn't deemed singular earlier,
             % it's highly unlikely. Could be numerical issue.
             % More robust check involves ranks as before, but costly.
             error('System potentially inconsistent or matrix ill-conditioned during forward substitution.');
         else
             rethrow(ME);
         end
    end
    
    % Solve Ux = y using backward substitution
    try
        x = linearSolvers.backwardSubstitution(U, y);
    catch ME
         if contains(ME.message, 'singular')
             % Check if U reveals inconsistency not caught earlier
             % (e.g., zero pivot with non-zero right side in last effective row)
             % A more rigorous check would look at rank(U) vs rank([U,y])
              error('System potentially inconsistent or matrix ill-conditioned during backward substitution.');
         else
             rethrow(ME);
         end
    end

end
