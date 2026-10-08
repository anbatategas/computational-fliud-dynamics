# Problem Set 5: Parabolic PDEs & The 1D Heat Equation

This final module focuses on transient (time-dependent) Partial Differential Equations. Specifically, it models the 1D diffusion/heat equation to track the spatio-temporal evolution of temperature across a conductive domain.

## The Mathematical Framework

The governing 1D parabolic PDE for diffusion is given by:

$$\frac{\partial T}{\partial t} = D \frac{\partial^2 T}{\partial x^2}$$

To solve this numerically over the domain $0 < x < L$ and $0 < t < T_{max}$, the spatial derivative is discretized using a second-order central difference. The time derivative is handled using three distinct time-marching schemes to evaluate their stability and accuracy.

## Implemented Solvers

All source codes are located in the `code/` directory. Because the implicit methods generate tridiagonal algebraic systems at every time step, the highly optimized `ThomasAlgorithmTriInput.m` direct solver is utilized to advance the solution through time efficiently.

### 1. Explicit Scheme (FTCS)
* **`Ex5_1_1_Explicit_Script.m`**: Implements Forward Time, Central Space. This method is computationally cheap per step but **conditionally stable**. It strictly requires the Fourier number (stability parameter $\lambda$) to remain below $0.5$ to prevent catastrophic numerical divergence.

### 2. Simple Implicit Scheme (BTCS)
* **`Ex5_1_2_Implicit_Script.m`**: Implements Backward Time, Central Space. This method is **unconditionally stable**, allowing for much larger time steps. It approaches the steady-state solution monotonically without any non-physical oscillations.

### 3. Crank-Nicolson Scheme
* **`Ex5_1_3_CrankNicolson_Script.m`**: A semi-implicit method that averages the spatial derivatives at the current and next time steps. It achieves $\mathcal{O}(\Delta t^2)$ accuracy in time (compared to the $\mathcal{O}(\Delta t)$ accuracy of the first two methods). While unconditionally stable, the numerical report demonstrates its tendency to produce slight numerical oscillations near sharp initial discontinuities before settling.

## Spatio-Temporal Evolution

The 3D surface plots below visualize the temperature decay over time across the 1D spatial domain, demonstrating the numerical dissipation of the initial conditions until a steady-state thermal equilibrium is reached.

<div style="display: flex; justify-content: center; align-items: center; gap: 10px;">
  <img src="images/Ex5_1_1_Fig1A.jpg" alt="Explicit Scheme 3D Evolution" width="48%">
  <img src="images/Ex5_1_3_Fig3A.jpg" alt="Crank-Nicolson Scheme 3D Evolution" width="48%">
</div>

*(Detailed 2D contour maps and constant-time profiles are available in the `images/` directory).*

## Documentation
For the complete stability condition derivations, truncation error proofs, and an in-depth analysis of the oscillatory nature of the Crank-Nicolson scheme, please refer to `CFD_PS5_Solutions.pdf`.
