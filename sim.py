import math
import matplotlib.pyplot as plt
import os

# --- Constants ---
g = 9.81
air_viscosity = 1.81e-5
air_density = 1.225

radius = 0.02135
mass = 0.04593
has_dimples = True

dt = 0.005

# --- Helper Functions ---

def calc_forces(vx, vy, vz, wx, wy, wz, omega):
    speed = math.sqrt(vx**2 + vy**2 + vz**2)
    if speed < 1e-6:
        speed = 1e-6

    Re = (air_density * speed * radius * 2) / air_viscosity
    S = omega * radius / speed

    #DO NOT TOUCH
    if has_dimples:
        CL = 1.99 * S - 3.25 * S**2  # quadruple check these equations
        
        # Quick double check of these equations, should be good to go.
        # Re at supercritical values follows linear fits that are based on spin rate
        # Re at lower values follows first order polynomial fit based on reynoldsNumber
        # Re at subcritical values based on separate fit, both from same literature
        if Re > 170000:
            C = 0.2 * S + 0.29
        elif Re > 150000:
             C = 0.2 * S + 0.28
        elif Re > 80000:
            C = 1.91e-11 * Re**2 - 5.40e-6 * Re + 0.56
        else:
            C = 1.29e-10 * Re * Re - 2.59e-5 * Re + 1.5
    else:
        CL = -0.38 * S**2 + 0.43 * S
        C = 0.5

    area = math.pi * radius**2
    liftMag = 0.5 * air_density * area * speed**2 * CL
    dragMag = 0.5 * air_density * area * C * speed**2

    crossX = vy * wz - vz * wy
    crossY = vz * wx - vx * wz
    crossZ = vx * wy - vy * wx
    crossMag = math.sqrt(crossX**2 + crossY**2 + crossZ**2) or 1e-9

    magnusForceX = liftMag * (crossX / crossMag)
    magnusForceY = liftMag * (crossY / crossMag)
    magnusForceZ = liftMag * (crossZ / crossMag)

    dragForceX = -dragMag * (vx / speed)
    dragForceY = -dragMag * (vy / speed)
    dragForceZ = -dragMag * (vz / speed)

    ax = (magnusForceX + dragForceX) / mass
    ay = (magnusForceY + dragForceY) / mass
    az = -g + (magnusForceZ + dragForceZ) / mass

    return ax, ay, az


def rk4_step(x, y, z, vx, vy, vz, wx, wy, wz, omega, dt):
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

    vx_new = vx + (k1vx + 2 * k2vx + 2 * k3vx + k4vx) / 6
    vy_new = vy + (k1vy + 2 * k2vy + 2 * k3vy + k4vy) / 6
    vz_new = vz + (k1vz + 2 * k2vz + 2 * k3vz + k4vz) / 6

    x_new = x + (k1x + 2 * k2x + 2 * k3x + k4x) / 6
    y_new = y + (k1y + 2 * k2y + 2 * k3y + k4y) / 6
    z_new = z + (k1z + 2 * k2z + 2 * k3z + k4z) / 6

    return x_new, y_new, z_new, vx_new, vy_new, vz_new


def simulate_shot(velocity, angle_deg, spin_rpm, shot_id=None):
    phirad = math.radians(angle_deg)
    vx = velocity * math.cos(phirad)
    vy = 0.0
    vz = velocity * math.sin(phirad)

    omega = spin_rpm * 2 * math.pi / 60
    wx, wy, wz = 0.0, omega, 0.0

    x, y, z = 0.0, 0.0, 0.0
    t = 0.0
    max_height = 0
    total_dist = 0

    # --- Store trajectory for plotting ---
    x_values, z_values = [x], [z]

    while z >= 0:
        x, y, z, vx, vy, vz = rk4_step(x, y, z, vx, vy, vz, wx, wy, wz, omega, dt)
        omega *= math.exp(-dt / 25)  # spin decay
        wy = omega
        t += dt

        #Convert m to yards and append data to be exported and analyzed
        x_values.append(x * 1.0936)
        z_values.append(z * 1.0936)

        if z > max_height:
            max_height = z
        if z <= 0 and x > 5:
            total_dist = x * 1.0936  # convert m to yards
            break

    # --- Plot and save trajectory ---
    plot_and_save_trajectory(x_values, z_values, shot_id=shot_id or "unknown")

    return total_dist, max_height, t


# --- Plotting Function ---
def plot_and_save_trajectory(x_values, z_values, 
                             shot_id=1,
                             title="Golf Ball Trajectory", 
                             x_label="Horizontal Distance (yards)", 
                             z_label="Height (yards)", 
                             color="blue",
                             save_dir="plots",
                             show_plot=False):
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
    #print(f"✅ Saved trajectory plot: {filepath}")


# --- Data --   -
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
    (63.7, 12.4, 2755, 234),
    (63.43, 12.3, 2778, 233),
    (65.04, 10, 2643, 231),
    (61.20, 13.4, 2813, 223),
    (58.47, 12, 2811, 209),
    (64.42, 6.1, 2686, 206),
    (57.89, 12.1, 2810, 206),
    (61.83, 7.2, 2817, 204),
    (63.79, 5.8, 2745, 204),
    (57.36, 7.1, 2808, 182),
    (53.60, 10.7, 2704, 179),
    (50.96, 11.6, 2580, 167),
    #Average of all shots
    (63.35, 10.88, 2679.62, 225.62)
]

# --- Simulation Loop ---
print(f"{'Velocity':>8} {'Angle':>6} {'Spin':>6} {'Sim Carry':>10} {'Expected':>10} {'Abs Error':>10} {'% Error'}")
for i, (velocity, angle, spin, expected) in enumerate(data, start=1):
    sim_carry, max_h, flight_time = simulate_shot(velocity, angle, spin, shot_id=i)
    error = sim_carry - expected
    errorPercent = 100 * error / expected
    print(f"{velocity:8.2f} {angle:6.2f} {spin:6} {sim_carry:10.2f} {expected:10} {error:10.2f} {errorPercent:10.2f}")
