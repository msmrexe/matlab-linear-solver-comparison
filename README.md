# MATLAB Linear Solver Comparison

This project implements and compares various direct and iterative methods for solving systems of linear equations ($Ax=b$) in MATLAB. It was developed for a Matrix Computations course to analyze the theoretical vs. practical performance of these fundamental algorithms using "from scratch" implementations where feasible.

The tool tests the following solvers:
* **Direct Methods (Implemented from Scratch):**
    * Gaussian Elimination (with partial pivoting & back-substitution)
    * Gauss-Jordan Elimination (with partial pivoting)
    * LU Decomposition (LUP Factorization, Forward & Backward Substitution)
* **Iterative Methods (Vectorized using Matrix Forms):**
    * Jacobi Method
    * Gauss-Seidel Method
    * Successive Over-Relaxation (SOR) (fixed omega=1.2)
* **Baseline:**
    * MATLAB's Backslash Operator (`A \ b`) (highly optimized direct solver)

It runs these solvers on matrices of varying sizes and types (random, diagonally dominant, ill-conditioned Hilbert), measures execution time and iteration counts (for iterative methods), and generates performance comparison plots.

## Features

* **Multiple Solvers:** Implements 6 core algorithms mostly from scratch (Jacobi, GS, SOR use `inv`; Convergence check uses `eig`) and compares against MATLAB's `\`.
* **"From Scratch" Direct Methods:** Gaussian, Gauss-Jordan, LUP Factorization, Forward/Backward Substitution are implemented using explicit loops and basic matrix indexing.
* **Vectorized Iterative Methods:** Jacobi, Gauss-Seidel, and SOR use efficient MATLAB matrix operations based on their theoretical matrix forms.
* **Robustness:** Includes partial pivoting for direct methods and convergence checks (spectral radius via `eig`) for iterative methods.
* **Diverse Testing:** Uses different matrix types (random, diagonally dominant, Hilbert) to show how matrix properties affect performance.
* **Performance Metrics:** Measures **Execution Time** and **Number of Iterations** (for iterative methods).
* **Clear Visualizations:** Generates `loglog` plots for time complexity (ideal for visualizing polynomial orders) and `semilogy` plots for iteration counts.
* **Modular Package:** Code is organized into a `+linearSolvers` package and utility functions.

## Project Structure

```
matlab-linear-solver-comparison/
├── .gitignore
├── README.md                         # This documentation
├── main.m                            # --- Run this script ---
├── analyzeSolvers.m                  # Script to run the analysis & save results
├── plotResults.m                     # Script to generate plots from results
├── +linearSolvers/                   # Package containing solver implementations
│   ├── gaussianElimination.m         # From-scratch Guassian Elimination
│   ├── gaussJordanElimination.m      # From-scratch Guass-Jordan Elimination
│   ├── luDecomposition.m             # Orchestrates LUP solving
│   ├── lupFactorization.m            # From-scratch LUP factorization
│   ├── forwardSubstitution.m         # From-scratch forward substitution
│   ├── backwardSubstitution.m        # From-scratch backward substitution
│   ├── jacobi.m                      # Vectorized Jacobi
│   ├── gaussSeidel.m                 # Vectorized Gauss-Seidel
│   ├── sor.m                         # Vectorized SOR
│   └── checkConvergence.m            # Utility for spectral radius check (uses eig)
└── utils/                            # Utility functions
  └── generateMatrices.m              # Function to create test matrices
```

## How the Algorithms Work

### Direct Methods
These methods aim to find the exact solution in a finite number of steps (ignoring floating-point errors). All direct methods here are implemented from scratch, including pivoting and substitution steps.

1.  **Gaussian Elimination:** Transforms the augmented matrix `[A|b]` into upper triangular form `[U|c]` using elementary row operations (forward elimination), then solves $Ux=c$ using back-substitution. Includes partial pivoting. Complexity: $O(n^3)$.
2.  **Gauss-Jordan Elimination:** Similar to Gaussian elimination, but continues row operations to transform `A` into the identity matrix `I`. The augmented matrix becomes `[I|x]`, directly revealing the solution $x$. Also uses pivoting. Complexity: $O(n^3)$.
3.  **LU Decomposition:** Factorizes matrix $A$ such that $PA = LU$, where $P$ is a permutation matrix (from pivoting), $L$ is lower triangular, and $U$ is upper triangular (using `lupFactorization.m`). Then solves $Ax=b$ by solving two simpler triangular systems: $Ly = Pb$ (`forwardSubstitution.m`) and $Ux=y$ (`backwardSubstitution.m`). Complexity: $O(n^3)$.

### Iterative Methods
These methods start with an initial guess $x^{(0)}$ and generate a sequence of approximations $x^{(1)}, x^{(2)}, \dots$ that hopefully converge to the true solution. They use vectorized matrix operations for efficiency. Convergence is checked using the spectral radius of the iteration matrix (calculated using MATLAB's `eig`).

1.  **Jacobi Method:** Matrix form: $x^{(k+1)} = D^{-1}(b - (L+U)x^{(k)})$.
2.  **Gauss-Seidel Method:** Matrix form: $x^{(k+1)} = (D+L)^{-1}(b - Ux^{(k)})$. (Uses MATLAB's `inv` for $(D+L)^{-1}$).
3.  **Successive Over-Relaxation (SOR):** Matrix form: $x^{(k+1)} = (D+\omega L)^{-1}(\omega b - (\omega U + (\omega-1)D)x^{(k)})$. (Uses MATLAB's `inv` for $(D+\omega L)^{-1}$).

**(Note:** While Jacobi is fully "from scratch" using basic operations, Gauss-Seidel and SOR use MATLAB's `inv` function for clarity and reasonable performance in the vectorized form. Implementing efficient iterative solvers for the triangular systems within GS/SOR is possible but adds significant complexity.)*

## How to Run

1.  **Open MATLAB.**
2.  **Navigate** to the project's root directory.
3.  **Run the main script** from the MATLAB Command Window:
    ```matlab
    >> main
    ```
4.  **Configuration:** Adjust parameters like `matrix_sizes`, `matrix_types`, `methods_to_test`, `tolerance`, and `max_iterations` directly within `main.m`.

The script will run the analysis, save results to `results/`, and generate plots in `plots/`.

## Performance Analysis

### Time vs. Size (LogLog Plot)

* **Direct Methods (Gaussian, Gauss-Jordan, LU, Backslash):** All these show lines roughly parallel to a slope of 3, indicating $O(n^3)$ complexity. MATLAB's `\` is the fastest. The from-scratch LU is slightly slower than `\` but competitive. Gaussian and Gauss-Jordan might are slightly slower still.
* **Iterative Methods (Jacobi, Gauss-Seidel, SOR):**
    * **For diagonally dominant matrices:** These to be faster than direct methods for larger $n$. Their lines have a slope closer to 2 (indicating roughly $O(n^2)$ work per iteration, with iteration count growing slowly). SOR/Gauss-Seidel generally outperform Jacobi.
    * **For random/Hilbert matrices:** These converge very slowly or fail (flagged in results). Their times are high, potentially hitting the iteration limit. Direct methods are more reliable here.

### Iterations vs. Size (SemiLogY Plot)

* **For diagonally dominant matrices:** Iteration counts remain low or grow very slowly.
* **For other matrices (random, Hilbert):** Iteration counts increase sharply with $n$, demonstrating the impact of poor conditioning on convergence.

---

## Author

Feel free to connect or reach out if you have any questions!

* **Maryam Rezaee**
* **GitHub:** [@msmrexe](https://github.com/msmrexe)
* **Email:** [ms.maryamrezaee@gmail.com](mailto:ms.maryamrezaee@gmail.com)

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for full details.
