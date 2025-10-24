% utils/generateMatrices.m
function A = generateMatrices(n, type)
% Generates different types of n x n matrices for testing.

    switch lower(type)
        case 'random'
            A = rand(n, n);
            % Ensure non-singularity (high probability with rand)
            if abs(det(A)) < 1e-10
                 A = generateMatrices(n, type); % Regenerate if singular
            end

        case 'diagdominant'
            A = rand(n, n);
            row_sums = sum(abs(A), 2) - abs(diag(A)); % Sum of off-diagonals
            % Add a value larger than row sums to the diagonal
            A = A + diag(row_sums + rand(n, 1) * n); % Make strictly diagonally dominant

        case 'hilbert'
            A = hilb(n); % MATLAB's built-in Hilbert matrix

        otherwise
            error('Unknown matrix type specified: %s', type);
    end
end
