# Computational Fluid Dynamics (CFD)

This repository contains MATLAB implementations of numerical methods, finite difference schemes, and computational solvers developed for the *Computational Fluid Dynamics* coursework at the National Technical University of Athens (NTUA). 

The repository is structured sequentially, progressing from fundamental linear algebra solvers and ODE integrators up to complex 2D spatial discretizations and transient partial differential equation (PDE) simulations.

## Repository Structure

Each module below has its own dedicated directory containing a specific `README.md` (detailing the mathematical framework), the MATLAB source code, execution scripts, and a comprehensive analytical report (PDF).

* **`Problem_Set_1/` : Numerical Solvers for Algebraic Systems**
  * Direct and iterative solvers for linear systems ($Ax=b$) and non-linear roots ($F(x)=0$). 
  * Features the Thomas Algorithm (TDMA) for tridiagonal matrices, alongside Jacobi, Gauss-Seidel, Successive Over-Relaxation (SOR), and Newton-Raphson methods.

* **`Problem_Set_2/` : Time Integration Methods for ODEs**
  * Explores stability and convergence of numerical integration schemes for Initial Value Problems. 
  * Compares 1st-order schemes (Explicit/Implicit Euler) with higher-order multi-stage schemes (Runge-Kutta 2nd, 3rd, and 4th order).

* **`Problem_Set_3/` : Boundary Value Problems & Poiseuille Flow**
  * Solves the 1D momentum equation (cylindrical coordinates) for steady, fully developed laminar flow in a pipe.
  * Discretizes the domain using second-order central finite differences and utilizes a highly memory-optimized Thomas Algorithm solver for the resulting tridiagonal system. Includes rigorous grid independence studies.

* **`Problem_Set_4/` : 2D Boundary Value Problems & Vectorized Solvers**
  * Expands to 2D elliptic PDEs by solving the Poisson equation for fully developed flow in a square duct using a 5-point central difference stencil.
  * **Highlight:** Features explicitly *vectorized* Jacobi and Gauss-Seidel solvers, drastically reducing computational overhead in MATLAB for dense ($101 \times 101$) grid iterations.

* **`Problem_Set_5/` : Parabolic PDEs & The 1D Heat Equation**
  * Models transient, time-dependent PDEs by tracking the spatio-temporal evolution of a 1D diffusion/heat equation.
  * Compares the conditional stability of the Explicit scheme (FTCS), the monotonic convergence of the Implicit scheme (BTCS), and the $\mathcal{O}(\Delta t^2)$ accuracy of the semi-implicit Crank-Nicolson scheme.

## Key Computational Skills Demonstrated
* **Finite Difference Methods:** 1D and 2D spatial discretization, boundary condition enforcement (including L'Hôpital's rule for singularities).
* **Algorithm Optimization:** Custom memory-efficient direct solvers (Tridiagonal Matrix Algorithms) and array vectorization for iterative methods.
* **Numerical Analysis:** Empirical verification of truncation errors, stability parameter ($\lambda$) constraints, and Chirikov resonance overlap criteria (via related Hamiltonian projects).

## Usage
All scripts are written in standard MATLAB and require no external toolboxes. Navigate to any `Problem_Set` directory and run the associated `Ex_#_Script.m` to generate the computational grids, solve the systems, and output the visualization plots.
