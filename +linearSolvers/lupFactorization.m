% +linearSolvers/lupFactorization.m
function [L, U, P] = lupFactorization(A)
% Performs LU decomposition with partial pivoting (LUP) from scratch.
% Factorizes PA = LU.

    [n, ~] = size(A);
    
    % Start with P as identity, U as a copy of A, L as identity
    P = eye(n);
    U = A; % U will be modified in place
    L = eye(n);

    for k = 1:n-1
        % --- Pivoting ---
        % Find the row with the largest element in column k, below or at row k
        [~, pivot_row_rel] = max(abs(U(k:n, k)));
        pivot_row_abs = pivot_row_rel + k - 1; % Absolute row index

        % Check for singularity
        if abs(U(pivot_row_abs, k)) < eps
             error('Matrix is singular or nearly singular (zero pivot found).');
        end

        % Swap rows in U, P, and the *part* of L already computed
        if pivot_row_abs ~= k
            U([k, pivot_row_abs], k:n) = U([pivot_row_abs, k], k:n); % Swap relevant part of U
            P([k, pivot_row_abs], :) = P([pivot_row_abs, k], :);     % Swap rows in P
            if k > 1
                L([k, pivot_row_abs], 1:k-1) = L([pivot_row_abs, k], 1:k-1); % Swap computed part of L
            end
        end

        % --- Elimination ---
        % Calculate multipliers and update L and U
        for i = k+1:n
            multiplier = U(i, k) / U(k, k);
            L(i, k) = multiplier; % Store multiplier in L
            % Update row i of U (perform elimination)
            U(i, k:n) = U(i, k:n) - multiplier * U(k, k:n);
            U(i, k) = 0; % Ensure it's exactly zero
        end
    end
     % Check the last pivot element after the loop
     if abs(U(n, n)) < eps
          error('Matrix is singular or nearly singular (zero pivot found).');
     end
end
