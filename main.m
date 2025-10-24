% main.m
% Main script to run the linear solver comparison analysis.

clear; clc; close all;
addpath(genpath(pwd)); % Add solvers and utils to path

% --- Configuration ---
matrix_sizes = [10, 20, 50, 100, 150]; % Sizes 'n' to test
matrix_types = {'random', 'diagdominant', 'hilbert'}; % Types of matrices
methods_to_test = { % Solvers to include in the comparison
    'GaussianElimination', ...
    'GaussJordanElimination', ...
    'LUDecomposition', ...
    'Jacobi', ...
    'GaussSeidel', ...
    'SOR', ...
    'MATLAB Backslash' ...
};
tolerance = 1e-6;       % Convergence tolerance for iterative methods
max_iterations = 2000;   % Max iterations for iterative methods

results_dir = 'results'; % Directory to save raw data
plots_dir = 'plots';     % Directory to save plots

results_filename = fullfile(results_dir, 'solver_results.mat');
% ---------------------

% --- Create output directories ---
if ~exist(results_dir, 'dir'); mkdir(results_dir); end
if ~exist(plots_dir, 'dir'); mkdir(plots_dir); end

% --- Run Analysis ---
fprintf('Starting Linear Solver Analysis...\n');
results_table = analyzeSolvers(matrix_sizes, matrix_types, methods_to_test, tolerance, max_iterations, results_filename);

% Display results table (optional)
disp('Analysis Complete. Results Summary:');
disp(results_table);

% --- Generate Plots ---
if ~isempty(results_table)
    fprintf('\nGenerating Plots...\n');
    plotResults(results_table, plots_dir);
else
    fprintf('\nNo results generated to plot.\n');
end

fprintf('\nAnalysis finished.\n');
