# Problem Set 1: Numerical Solvers for Algebraic Systems

This module focuses on the foundational numerical methods used to solve systems of equations that frequently arise when discretizing partial differential equations (PDEs) in fluid mechanics.

## Implemented Methods

### 1. Thomas Algorithm (TDMA)
An efficient direct solver for tridiagonal systems of linear equations, heavily used in 1D implicit finite difference schemes.
* `ThomasAlgorithm.m`: Core algorithm implementation.
* `ThomasAlgorithmEx1_1_Script.m`: Solution script for Exercise 1.1.

### 2. Iterative Solvers for Linear Systems ($Ax=b$)
Algorithms to solve diagonally dominant square matrices. The implementations utilize weighted infinity norms as termination criteria.
* `JacobiMethod.m`: Jacobi iteration.
* `GaussSeidelMethod.m`: Gauss-Seidel iteration.
* `SORMethod.m`: Successive Over-Relaxation (SOR) with dynamic relaxation parameter $\omega$.
* Execution scripts: `JacobiMethodEx1_2_Script.m`, `GaussSeidelMethodEx1_2_Script.m`, `SORMethodEx1_2_Script.m`.

### 3. Newton-Raphson Method for Non-Linear Systems ($F(x)=0$)
Root-finding algorithms for non-linear multi-variable systems. 
* `NRMethod.m`: Computes the Jacobian matrix numerically via finite differences.
* `NRMethod_Jin.m`: Accepts an analytically derived Jacobian matrix as a function handle.
* Execution scripts: `NRMethodEx1_3_Script.m`, `NRMethod_JinEx1_3_Script.m`.

## Documentation
The rigorous mathematical formulation, matrix stability checks, and convergence analyses can be found in `CFD_PS1_Solutions.pdf`.
