"""
Golf Ball Flight Simulation
---------------------------
This Python module implements a reproducible, physics-based model of golf ball
flight using empirically derived aerodynamic relations. The simulation resolves
translational motion under gravity, drag, and lift (Magnus) forces, while
accounting for Reynolds number and spin ratio dependencies in both the dimpled
and smooth-sphere regimes.

The model uses:
    - Fourth-order Runge–Kutta (RK4) numerical integration
    - Empirical C_L(S) and C_D(Re, S) relations
    - Exponential spin decay over time

References:
    • Smits & Smith (1994) — aerodynamic data for dimpled spheres
    • Muto et al. (2012) — smooth-sphere wind-tunnel correlations
    • Davies (1949) — classical smooth-sphere flow characterization
"""

import math
import matplotlib.pyplot as plt
import os

# ----------------------------------------------------------------------
# Global Constants
# ----------------------------------------------------------------------

g = 9.81                     # gravitational acceleration [m/s²]
air_viscosity = 1.81e-5      # dynamic viscosity of air [Pa·s]
air_density = 1.225          # air density at sea level [kg/m³]

radius = 0.02135             # golf ball radius [m]
mass = 0.04593               # golf ball mass [kg]
has_dimples = True           # surface condition toggle (True = dimpled)
dt = 0.005                   # integration time step [s]


# ----------------------------------------------------------------------
# Aerodynamic Force Computation
# ----------------------------------------------------------------------

def calc_forces(vx, vy, vz, wx, wy, wz, omega):
    """
    Compute instantaneous aerodynamic accelerations.

    Parameters
    ----------
    vx, vy, vz : float
        Cartesian velocity components [m/s].
    wx, wy, wz : float
        Angular velocity components [rad/s].
    omega : float
        Total angular speed [rad/s].

    Returns
    -------
    ax, ay, az : float
        Accelerations [m/s²] including drag, lift (Magnus), and gravity.
    """
    speed = math.sqrt(vx**2 + vy**2 + vz**2)
    if speed < 1e-6:  # avoid division by zero
        speed = 1e-6

    # Dimensionless aerodynamic parameters
    Re = (air_density * speed * radius * 2) / air_viscosity  # Reynolds number
    S = omega * radius / speed                               # spin ratio

    # ------------------------------------------------------------------
    # Aerodynamic Coefficient Models
    # ------------------------------------------------------------------
    # This section implements piecewise C_D(Re, S) and C_L(S) relations
    # for both dimpled and smooth golf balls. Each model reproduces
    # empirical data from wind-tunnel literature and driver-class launch
    # conditions.

    if has_dimples:
        # ---- Dimpled Ball: Lift Model ----
        # Empirical quadratic fit from dimpled-ball data showing lift rise
        # and eventual decay due to partial flow separation at high spin.
        CL = 1.99 * S - 3.25 * S**2

        # ---- Dimpled Ball: Drag Model ----
        # Piecewise C_D(Re, S) representation spanning laminar,
        # transitional, and supercritical regimes.
        if Re > 170000:
            C = 0.2 * S + 0.29
        elif Re > 150000:
            C = 0.2 * S + 0.28
        elif Re > 80000:
            C = 1.91e-11 * Re**2 - 5.40e-6 * Re + 0.56
        else:
            C = 1.29e-10 * Re**2 - 2.59e-5 * Re + 1.5

    else:
        # ---- Smooth Ball: Lift Model ----
        # Smooth spheres produce much weaker lift; two polynomial regimes
        # (above/below Re = 2×10⁵) based on Muto et al. (2012).
        Re5 = Re / 1e5
        a = 8.838384
        if Re > 200000:
            b = -3.017677   # High-Re: diminished slope
        else:
            b = -2.0        # Low-Re: enhanced slope
        CL = a * S**2 + b * S

        # ---- Smooth Ball: Drag Model ----
        # Polynomial from Muto et al. (2012) for smooth spheres
        # in 0.8–2.0×10⁵ Re range, capped at C=0.48 for stability.
        a = -0.115786
        b = 0.013145
        c = 0.740801
        C = min(a * Re5**2 + b * Re5 + c, 0.48)

    # ------------------------------------------------------------------
    # Force Calculations
    # ------------------------------------------------------------------
    area = math.pi * radius**2

    # Lift and drag magnitudes (aerodynamic forces)
    liftMag = 0.5 * air_density * area * speed**2 * CL
    dragMag = 0.5 * air_density * area * C * speed**2

    # Direction of Magnus (lift) force via cross product of v × ω
    crossX = vy * wz - vz * wy
    crossY = vz * wx - vx * wz
    crossZ = vx * wy - vy * wx
    crossMag = max(math.sqrt(crossX**2 + crossY**2 + crossZ**2), 1e-9)

    magnusForceX = liftMag * (crossX / crossMag)
    magnusForceY = liftMag * (crossY / crossMag)
    magnusForceZ = liftMag * (crossZ / crossMag)

    # Direction of drag force (opposite velocity)
    dragForceX = -dragMag * (vx / speed)
    dragForceY = -dragMag * (vy / speed)
    dragForceZ = -dragMag * (vz / speed)

    # Resultant accelerations
    ax = (magnusForceX + dragForceX) / mass
    ay = (magnusForceY + dragForceY) / mass
    az = -g + (magnusForceZ + dragForceZ) / mass

    return ax, ay, az


# ----------------------------------------------------------------------
# Numerical Integration: Fourth-Order Runge–Kutta Scheme
# ----------------------------------------------------------------------

def rk4_step(x, y, z, vx, vy, vz, wx, wy, wz, omega, dt):
    """
    Advance the system state by one time step using a 4th-order
    Runge–Kutta integrator.

    Returns updated position and velocity components after Δt.
    """
    ax1, ay1, az1 = calc_forces(vx, vy, vz, wx, wy, wz, omega)
    k1vx, k1vy, k1vz = ax1 * dt, ay1 * dt, az1 * dt
    k1x, k1y, k1z = vx * dt, vy * dt, vz * dt

    ax2, ay2, az2 = calc_forces(vx + 0.5 * k1vx, vy + 0.5 * k1vy, vz + 0.5 * k1vz, wx, wy, wz, omega)
    k2vx, k2vy, k2vz = ax2 * dt, ay2 * dt, az2 * dt
    k2x, k2y, k2z = (vx + 0.5 * k1vx) * dt, (vy + 0.5 * k1vy) * dt, (vz + 0.5 * k1vz) * dt

    ax3, ay3, az3 = calc_forces(vx + 0.5 * k2vx, vy + 0.5 * k2vy, vz + 0.5 * k2vz, wx, wy, wz, omega)
    k3vx, k3vy, k3vz = ax3 * dt, ay3 * dt, az3 * dt
    k3x, k3y, k3z = (vx + 0.5 * k2vx) * dt, (vy + 0.5 * k2vy) * dt, (vz + 0.5 * k2vz) * dt

    ax4, ay4, az4 = calc_forces(vx + k3vx, vy + k3vy, vz + k3vz, wx, wy, wz, omega)
    k4vx, k4vy, k4vz = ax4 * dt, ay4 * dt, az4 * dt
    k4x, k4y, k4z = (vx + k3vx) * dt, (vy + k3vy) * dt, (vz + k3vz) * dt

    # Weighted average (RK4)
    vx_new = vx + (k1vx + 2 * k2vx + 2 * k3vx + k4vx) / 6
    vy_new = vy + (k1vy + 2 * k2vy + 2 * k3vy + k4vy) / 6
    vz_new = vz + (k1vz + 2 * k2vz + 2 * k3vz + k4vz) / 6

    x_new = x + (k1x + 2 * k2x + 2 * k3x + k4x) / 6
    y_new = y + (k1y + 2 * k2y + 2 * k3y + k4y) / 6
    z_new = z + (k1z + 2 * k2z + 2 * k3z + k4z) / 6

    return x_new, y_new, z_new, vx_new, vy_new, vz_new


# ----------------------------------------------------------------------
# Shot Simulation (Trajectory Integration)
# ----------------------------------------------------------------------

def simulate_shot(velocity, angle_deg, spin_rpm, shot_id=None):
    """
    Simulate a single golf shot and compute trajectory metrics.

    Parameters
    ----------
    velocity : float
        Launch velocity [m/s]
    angle_deg : float
        Launch angle [degrees]
    spin_rpm : float
        Initial backspin rate [RPM]

    Returns
    -------
    total_dist : float
        Horizontal carry distance [yards]
    max_height : float
        Maximum height [m]
    t : float
        Flight duration [s]
    """
    phirad = math.radians(angle_deg)
    vx, vy, vz = velocity * math.cos(phirad), 0.0, velocity * math.sin(phirad)

    omega = spin_rpm * 2 * math.pi / 60  # convert RPM → rad/s
    wx, wy, wz = 0.0, omega, 0.0

    x = y = z = 0.0
    t = 0.0
    max_height = 0.0
    total_dist = 0.0

    x_values, z_values = [x], [z]

    while z >= 0:
        x, y, z, vx, vy, vz = rk4_step(x, y, z, vx, vy, vz, wx, wy, wz, omega, dt)
        omega *= math.exp(-dt / 25)  # spin decay (τ ≈ 25 s)
        wy = omega
        t += dt

        x_values.append(x * 1.0936)
        z_values.append(z * 1.0936)

        if z > max_height:
            max_height = z
        if z <= 0 and x > 5:
            total_dist = x * 1.0936
            break

    return total_dist, max_height, t

# ----------------------------------------------------------------------
# Plotting Utilities
# ----------------------------------------------------------------------

def plot_and_save_trajectory(x_values, z_values, 
                             shot_id=1,
                             title="Golf Ball Trajectory", 
                             x_label="Horizontal Distance (yards)", 
                             z_label="Height (yards)", 
                             color="blue",
                             save_dir="plots",
                             show_plot=False):
    """
    Plot and export a single shot trajectory.

    Parameters
    ----------
    x_values, z_values : list of float
        Horizontal and vertical positions (in yards) along the trajectory.
    shot_id : int
        Identifier for labeling the saved plot.
    title : str
        Plot title.
    x_label, z_label : str
        Axis labels.
    color : str
        Trajectory color.
    save_dir : str
        Directory in which to save PNG files.
    show_plot : bool
        If True, display the plot interactively; otherwise, save silently.

    Notes
    -----
    The figure is styled for clarity with gridlines and proportional aspect
    ratio, allowing consistent comparison between trajectories.
    """
    os.makedirs(save_dir, exist_ok=True)
    plt.figure(figsize=(8, 4))
    plt.plot(x_values, z_values, color=color, linewidth=2)
    plt.title(f"{title} (Shot {shot_id})")
    plt.xlabel(x_label)
    plt.ylabel(z_label)
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()

    filepath = os.path.join(save_dir, f"trajectory_shot{shot_id}.png")
    plt.savefig(filepath, dpi=150)

    if show_plot:
        plt.show()
    else:
        plt.close()


# ----------------------------------------------------------------------
# Automated Trajectory Generator
# ----------------------------------------------------------------------

def generate_trajectory(velocity, angle_deg, spin_rpm, 
                        has_dimples=False, shot_id=None, 
                        save_dir="plots", dt=0.005):
    """
    Generate, integrate, and save a full 2D trajectory image for a single shot.

    Parameters
    ----------
    velocity : float
        Launch velocity [m/s].
    angle_deg : float
        Launch angle [degrees].
    spin_rpm : float
        Initial backspin [rpm].
    has_dimples : bool
        Surface condition flag.
    shot_id : int
        Optional identifier for file naming.
    save_dir : str
        Directory in which to save trajectory plots.
    dt : float
        Integration time step [s].

    Returns
    -------
    total_dist : float
        Carry distance [yards].
    max_height : float
        Maximum trajectory height [yards].
    t : float
        Flight duration [s].

    Notes
    -----
    The function internally calls the RK4 integrator and applies exponential
    spin decay to simulate angular momentum loss during flight. The output plot
    mirrors standard launch-monitor trajectory displays.
    """
    phirad = math.radians(angle_deg)
    vx = velocity * math.cos(phirad)
    vy = 0.0
    vz = velocity * math.sin(phirad)

    omega = spin_rpm * 2 * math.pi / 60
    wx, wy, wz = 0.0, omega, 0.0

    x, y, z = 0.0, 0.0, 0.0
    t = 0.0
    max_height = 0.0
    total_dist = 0.0

    # Initialize storage for trajectory visualization
    x_values, z_values = [x], [z]

    while z >= 0:
        x, y, z, vx, vy, vz = rk4_step(x, y, z, vx, vy, vz, wx, wy, wz, omega, dt)
        omega *= math.exp(-dt / 25)  # exponential spin decay (τ ≈ 25 s)
        wy = omega
        t += dt

        # Record trajectory in yards
        x_values.append(x * 1.0936)
        z_values.append(z * 1.0936)

        if z > max_height:
            max_height = z
        if z <= 0 and x > 5:
            total_dist = x * 1.0936
            break

    # --- Export Trajectory Plot ---
    os.makedirs(save_dir, exist_ok=True)
    plt.figure(figsize=(8, 4))
    plt.plot(x_values, z_values, color="blue", linewidth=2)
    plt.title(f"Trajectory (v={velocity:.1f} m/s, θ={angle_deg:.1f}°, ω={spin_rpm} rpm)")
    plt.xlabel("Horizontal Distance (yards)")
    plt.ylabel("Height (yards)")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()

    filepath = os.path.join(save_dir, f"trajectory_{shot_id or 'auto'}.png")
    plt.savefig(filepath, dpi=150)
    plt.close()

    return total_dist, max_height * 1.0936, t


# ----------------------------------------------------------------------
# Experimental / Reference Data
# ----------------------------------------------------------------------
# Each tuple: (velocity [m/s], launch angle [°], spin [rpm], measured carry [yd])
# Data sources:
#   - Indoor simulator measurements (dimpled)
#   - Literature-referenced driver-class data (smooth & dimpled)
#   - 2018 Proceedings study; Muto et al. (2012); McIlroy reference shot
data = [
    (70.68, 12.4, 2507, 268),
    (69.69, 10.5, 2539, 257),
    (67.32, 14.1, 2603, 256),
    (66.92, 12.4, 2612, 250),
    (67.77, 11.1, 2592, 250),
    (67.28, 12.3, 2604, 247),
    (66.03, 12.3, 2625, 245),
    (66.25, 12.4, 2623, 245),
    (66.70, 10.3, 2617, 242),
    (63.70, 12.4, 2755, 234),
    (63.43, 12.3, 2778, 233),
    (65.04, 10.0, 2643, 231),
    (61.20, 13.4, 2813, 223),
    (58.47, 12.0, 2811, 209),
    (64.42, 6.1, 2686, 206),
    (57.89, 12.1, 2810, 206),
    (61.83, 7.2, 2817, 204),
    (63.79, 5.8, 2745, 204),
    (57.36, 7.1, 2808, 182),
    (53.60, 10.7, 2704, 179),
    (50.96, 11.6, 2580, 167),
    (71.50, 6.5, 3000, 285),  # 2018 Proceedings (dimpled)
    (82.00, 14.0, 2400, 109), # Estimated smooth-ball drive (R. McIlroy)
]


# ----------------------------------------------------------------------
# Batch Simulation Loop and Performance Summary
# ----------------------------------------------------------------------

print(f"{'Velocity':>8} {'Angle':>6} {'Spin':>6} {'Sim Carry':>10} "
      f"{'Expected':>10} {'Abs Error':>10} {'% Error'}")

for i, (velocity, angle, spin, expected) in enumerate(data, start=1):
    # Run simulation (RK4 integration + aerodynamic modeling)
    sim_carry, max_h, flight_time = generate_trajectory(velocity, angle, spin, shot_id=i)

    # Compare model output with experimental or literature distances
    error = sim_carry - expected
    errorPercent = 100 * error / expected

    # Formatted output for validation and reporting
    print(f"{velocity:8.2f} {angle:6.2f} {spin:6} "
          f"{sim_carry:10.2f} {expected:10} {error:10.2f} {errorPercent:10.2f}")
