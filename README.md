# Golf Ball Flight Simulation

### Overview
This repository contains a reproducible, physics-based simulation framework for modeling golf ball flight under varying aerodynamic conditions.  
The code implements a fourth-order Runge–Kutta (RK4) numerical integrator and empirically derived aerodynamic models that account for the effects of spin, Reynolds number, and surface texture (dimpled vs. smooth).

The project serves as both a **computational study** and a **validation tool**, comparing simulated trajectories against launch monitor data and published wind-tunnel measurements. The model captures the distinct aerodynamic behavior of dimpled and smooth golf balls, emphasizing their differing drag crises and lift generation characteristics.

---

### Features
- **Empirical Aerodynamics**
  - Lift and drag coefficients derived from **Smits & Smith (1994)** and **Muto et al. (2012)**.
  - Piecewise drag model for dimpled spheres, spanning laminar to supercritical regimes.
  - Reynolds-dependent lift model for smooth spheres with high- and low-Reynolds fits.

- **Numerical Integration**
  - Full 3D translational motion solved via a fourth-order **Runge–Kutta (RK4)** scheme.
  - Continuous recalculation of aerodynamic forces at each time step.
  - Exponential spin decay to represent angular momentum loss.

- **Validation**
  - Includes experimental and simulator carry distances for 21 shots.
  - Produces comparative metrics: simulated vs. measured carry distance, absolute error, and percent error.
  - Automatically exports trajectory plots for each shot in a clean, publication-ready format.

---

### Equations and Models

**Lift Coefficient (Dimpled Ball)**  
C_L = 1.99S − 3.25S²


**Lift Coefficient (Smooth Ball)**  
C_L = 8.833S² − 3.017S    (for Re > 2×10⁵)
C_L = 8.833S² − 2.000S    (for Re ≤ 2×10⁵)


**Drag Coefficient (Dimpled Ball)**  
C_D(Re, S) = {
    0.2·S + 0.29                         for Re > 170,000
    0.2·S + 0.28                         for 150,000 < Re ≤ 170,000
    1.91×10⁻¹¹·Re² − 5.40×10⁻⁶·Re + 0.56 for 80,000 < Re ≤ 150,000
    1.29×10⁻¹⁰·Re² − 2.59×10⁻⁵·Re + 1.50 for Re ≤ 80,000
}


**Drag Coefficient (Smooth Ball)**  
C_D = min( −0.1158(Re/10⁵)² + 0.0131(Re/10⁵) + 0.7408 , 0.48 )


---

### References
1. Smits, A. J., & Smith, D. R. (1994). *A new aerodynamic model of a golf ball in flight.*  
   *Science and Golf II: Proceedings of the World Scientific Congress of Golf*, E & F.N. Spon, London.

2. Muto, M., Nakayama, T., & Watanabe, K. (2012). *Experimental and numerical analysis of the aerodynamic characteristics of a smooth and dimpled sphere.*  
   *Procedia Engineering*, 34, 125–130.

3. Davies, C. N. (1949). *The aerodynamics of large spheres at low and moderate Reynolds numbers.*  
   *Proc. Phys. Soc. London, Series B*, 62(5), 259–270.

---

### Outputs
- Trajectory plots (`/plots/trajectory_shot#.png`)
- Tabulated carry-distance comparison (terminal output)
- Optional 3D surface plots for coefficient visualizations (e.g., \( C_L(S, Re) \))

Example console output:

