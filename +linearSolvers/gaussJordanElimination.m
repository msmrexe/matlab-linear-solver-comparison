% +linearSolvers/gaussJordanElimination.m
function x = gaussJordanElimination(A, b)
% Solves Ax = b using Gauss-Jordan Elimination with partial pivoting.

    [n, ~] = size(A);
    Ab = [A, b(:)]; % Augmented matrix

    % Forward Elimination (similar to Gaussian, but reduces above pivot too)
    for i = 1:n
        % Find pivot row (only search below and including current row)
        [~, maxRowIdx] = max(abs(Ab(i:n, i)));
        maxRowIdx = maxRowIdx + i - 1;

        % Swap rows
        if maxRowIdx ~= i
            Ab([i, maxRowIdx], :) = Ab([maxRowIdx, i], :);
        end

        % Check for singular matrix
        pivot = Ab(i, i);
        if abs(pivot) < eps
            % Check for inconsistency before declaring singular
             if any(abs(Ab(i:n, n+1)) > eps)
                 error('System has no solution (inconsistent).');
             else
                 % This implies infinite solutions if rank < n,
                 % but Gauss-Jordan aims for a unique solution.
                 error('Matrix is singular or nearly singular. Cannot find unique solution.');
             end
        end

        % Normalize the pivot row
        Ab(i, :) = Ab(i, :) / pivot;

        % Eliminate column i in ALL other rows (above and below)
        for j = [1:i-1, i+1:n] % Rows other than i
            factor = Ab(j, i); % Element to eliminate
            Ab(j, :) = Ab(j, :) - factor * Ab(i, :);
        end
    end

    % Solution is the last column
    x = Ab(:, n+1);
end
