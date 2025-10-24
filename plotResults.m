% plotResults.m
% Generates plots from the linear solver analysis results.

function plotResults(results_table, output_dir)

    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    methods = unique(results_table.Method);
    matrix_types = unique(results_table.MatrixType);
    markers = {'-o', '-s', '-^', '-d', '-v', '-p', '-h', '-x'}; % Marker styles

    % --- Plot 1: Time vs. Size (LogLog Scale) ---
    figure('Name', 'Time Complexity Analysis', 'NumberTitle', 'off');
    hold on;
    xlabel('Matrix Size (n)');
    ylabel('Execution Time (s)');
    title('Solver Time Complexity (LogLog Scale)');
    set(gca, 'XScale', 'log', 'YScale', 'log'); % Use loglog plot
    grid on;
    legend_entries_time = {};

    marker_idx = 1;
    for i = 1:length(methods)
        method_name = methods{i};
        for j = 1:length(matrix_types)
            m_type = matrix_types{j};
            
            % Filter data for this method and matrix type
            subset = results_table(strcmp(results_table.Method, method_name) & ...
                                   strcmp(results_table.MatrixType, m_type) & ...
                                   results_table.Converged == true, :); % Only plot converged/successful runs

            if ~isempty(subset)
                plot(subset.Size, subset.Time_s, markers{marker_idx}, 'LineWidth', 1.5, 'MarkerSize', 6);
                legend_entries_time{end+1} = sprintf('%s (%s)', method_name, m_type);
                 marker_idx = mod(marker_idx, length(markers)) + 1;
            end
        end
    end
    legend(legend_entries_time, 'Location', 'northwest');
    hold off;
    saveas(gcf, fullfile(output_dir, 'time_vs_size_loglog.png'));
    
    % --- Plot 2: Iterations vs. Size (SemiLogY Scale) ---
    iterative_methods_results = results_table(~isnan(results_table.Iterations), :);
    if ~isempty(iterative_methods_results)
        iterative_methods = unique(iterative_methods_results.Method);
        
        figure('Name', 'Iteration Analysis', 'NumberTitle', 'off');
        hold on;
        xlabel('Matrix Size (n)');
        ylabel('Number of Iterations');
        title('Convergence Rate (Iterations vs. Size)');
         set(gca, 'YScale', 'log'); % Use semilogy for iterations
        grid on;
        legend_entries_iter = {};
        
        marker_idx = 1;
         for i = 1:length(iterative_methods)
             method_name = iterative_methods{i};
              for j = 1:length(matrix_types)
                 m_type = matrix_types{j};
                 
                 subset = iterative_methods_results(strcmp(iterative_methods_results.Method, method_name) & ...
                                                     strcmp(iterative_methods_results.MatrixType, m_type) & ...
                                                     iterative_methods_results.Converged == true, :);
                 
                  if ~isempty(subset)
                     plot(subset.Size, subset.Iterations, markers{marker_idx}, 'LineWidth', 1.5, 'MarkerSize', 6);
                     legend_entries_iter{end+1} = sprintf('%s (%s)', method_name, m_type);
                      marker_idx = mod(marker_idx, length(markers)) + 1;
                  end
              end
         end
        legend(legend_entries_iter, 'Location', 'northwest');
        hold off;
        saveas(gcf, fullfile(output_dir, 'iterations_vs_size_semilogy.png'));
    else
        disp('No iteration data to plot.');
    end

    fprintf('Plots saved to directory: %s\n', output_dir);
end
