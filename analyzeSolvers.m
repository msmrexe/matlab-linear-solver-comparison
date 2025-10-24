% analyzeSolvers.m
% Runs timing and iteration analysis for linear solvers.

function results_table = analyzeSolvers(sizes, matrix_types, methods_to_test, tol, max_iter, results_file)
    % Add package path
    addpath(genpath(pwd));

    % Define which methods are iterative
    iterative_methods = {'Jacobi', 'GaussSeidel', 'SOR'};

    results = []; % Use a cell array to accumulate results

    for n = sizes
        fprintf('Analyzing matrix size n = %d...\n', n);
        for m_type = matrix_types
            fprintf('  Matrix Type: %s...\n', m_type);

            % Generate matrix and vector
            try
                A = utils.generateMatrices(n, m_type);
                b = rand(n, 1);
            catch ME
                fprintf('    Skipping matrix generation error: %s\n', ME.message);
                continue;
            end

            for method_name_cell = methods_to_test
                method_name = method_name_cell{:}; % Get string from cell
                fprintf('    Testing Method: %s...\n', method_name);

                time_taken = NaN;
                iterations = NaN; % Only for iterative methods
                converged = true; % Assume success unless failure occurs
                error_msg = '';

                try
                    % Check applicability/convergence for iterative methods
                    is_iterative = ismember(method_name, iterative_methods);
                    can_run = true;
                    if is_iterative
                        % Use spectral radius check for convergence guarantee
                        % Using omega=1.2 as a default for SOR check if needed
                        omega_check = 1.2;
                         if ~linearSolvers.checkConvergence(A, method_name, omega_check)
                             warning('Skipping %s for n=%d (%s matrix): Convergence condition (spectral radius >= 1) not met.', method_name, n, m_type);
                             can_run = false;
                             error_msg = 'Convergence Fail';
                         end
                    end

                    if can_run
                        tic; % Start timer
                        switch method_name
                            case 'GaussianElimination'
                                x = linearSolvers.gaussianElimination(A, b);
                            case 'GaussJordanElimination'
                                x = linearSolvers.gaussJordanElimination(A, b);
                            case 'LUDecomposition'
                                x = linearSolvers.luDecomposition(A, b);
                            case 'Jacobi'
                                [x, iterations, final_error] = linearSolvers.jacobi(A, b, [], tol, max_iter);
                                if final_error >= tol; converged = false; error_msg = sprintf('Did not converge (err=%e)', final_error); end
                            case 'GaussSeidel'
                                [x, iterations, final_error] = linearSolvers.gaussSeidel(A, b, [], tol, max_iter);
                                if final_error >= tol; converged = false; error_msg = sprintf('Did not converge (err=%e)', final_error); end
                            case 'SOR'
                                omega = 1.2; % Example omega, could be optimized
                                [x, iterations, final_error] = linearSolvers.sor(A, b, omega, [], tol, max_iter);
                                if final_error >= tol; converged = false; error_msg = sprintf('Did not converge (err=%e)', final_error); end
                            case 'MATLAB Backslash'
                                x = A \ b;
                        end
                        time_taken = toc; % Stop timer
                    end

                catch ME
                    warning('Error running %s for n=%d (%s matrix): %s', method_name, n, m_type, ME.message);
                    error_msg = ME.identifier; % Record error type
                    converged = false;
                end

                % Append result
                results = [results; {method_name, n, m_type, time_taken, iterations, converged, error_msg}];
            end % End methods loop
        end % End matrix types loop
    end % End sizes loop

    % Convert cell array to table
    results_table = cell2table(results, ...
        'VariableNames', {'Method', 'Size', 'MatrixType', 'Time_s', 'Iterations', 'Converged', 'Error'});

    % Save results
    if nargin >= 6 && ~isempty(results_file)
        try
            save(results_file, 'results_table');
            fprintf('\nResults saved to %s\n', results_file);
        catch ME
            fprintf('\nError saving results: %s\n', ME.message);
        end
    end
end
