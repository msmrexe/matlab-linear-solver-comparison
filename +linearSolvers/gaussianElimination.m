% +linearSolvers/gaussianElimination.m
function x = gaussianElimination(A, b)
% Solves Ax = b using Gaussian Elimination with partial pivoting
% followed by back-substitution.

    [n, ~] = size(A);
    Ab = [A, b(:)]; % Augmented matrix

    % Forward Elimination with Partial Pivoting
    for i = 1:n-1
        % Find pivot row
        [~, maxRowIdx] = max(abs(Ab(i:n, i)));
        maxRowIdx = maxRowIdx + i - 1; % Adjust index relative to full matrix

        % Swap rows if necessary
        if maxRowIdx ~= i
            Ab([i, maxRowIdx], :) = Ab([maxRowIdx, i], :);
        end

        % Check for singular matrix (zero pivot after pivoting)
        if abs(Ab(i, i)) < eps % Use epsilon for floating point comparison
            error('Matrix is singular or nearly singular. Cannot solve using Gaussian Elimination.');
        end

        % Eliminate column i in rows below i
        for j = i+1:n
            factor = Ab(j, i) / Ab(i, i);
            Ab(j, i:n+1) = Ab(j, i:n+1) - factor * Ab(i, i:n+1);
        end
    end

    % Check for inconsistency in the last row
    if abs(Ab(n, n)) < eps && abs(Ab(n, n+1)) > eps
         error('System has no solution (inconsistent).');
    elseif abs(Ab(n, n)) < eps && abs(Ab(n, n+1)) < eps
         error('System has infinitely many solutions.'); % Or handle appropriately
    end

    % Back Substitution (Pure MATLAB, no SymPy)
    x = zeros(n, 1);
    x(n) = Ab(n, n+1) / Ab(n, n);
    for i = n-1:-1:1
        sum_ax = Ab(i, i+1:n) * x(i+1:n);
        x(i) = (Ab(i, n+1) - sum_ax) / Ab(i, i);
    end

end
