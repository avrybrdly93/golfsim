import math
import js

# Get the JS canvas context
ctx = js.canvas.getContext("2d")

# ========== INITIAL SETUP ==========
firstRun = True
spinRate = 3000  # RPM
phideg = 15
initialVelocity = 70  # m/s
airDensity = 1.225  # kg/m^3
airViscosity = 1.81e-5
radius = 0.02135  # m (golf ball)
dt = 0.01  # timestep
hasDimples = True

# Gravity
g = 9.81

# Initial conditions
velocityY = 0
wx, wy, wz = 0, 0, 0
x, y, z = 0, 0, 1
t = 0
count = 0

# Lists for plotting
trajectory = []

# ========== FIRST PASS ==========
if firstRun:
    omega = spinRate * 2 * math.pi / 60.0
    print("This is the first time!")
    phirad = math.pi * phideg / 180
    velocity = initialVelocity
    velocityX = initialVelocity * math.cos(phirad)
    velocityZ = initialVelocity * math.sin(phirad)
    print(f"First evolution, spinrate, launch angle, and initial velocity: {omega:.2f}, {velocity:.2f}, {phideg}")
    firstRun = False

# ========== SIMULATION LOOP ==========
while z > 0:
    # Calculate current conditions
    speed = math.sqrt(velocityX**2 + velocityY**2 + velocityZ**2)
    safeV = max(speed, 0.1)
    reynoldsNumber = (airDensity * speed * radius * 2) / airViscosity
    S = (omega * radius) / safeV
    omega *= math.exp(-dt / 4.5)  # spin decay

    if hasDimples:
        CL = 1.99 * S - 3.25 * S**2
        if reynoldsNumber > 175000:
            C = min(1.91e-11 * reynoldsNumber**2 - 5.40e-6 * reynoldsNumber + 0.56, 0.21)
        elif reynoldsNumber > 80000:
            C = 1.91e-11 * reynoldsNumber**2 - 5.40e-6 * reynoldsNumber + 0.56
        else:
            C = 0.5
    else:
        C = 0.5
        CL = -0.38 * S**2 + 0.43 * S

    area = math.pi * radius**2

    # Forces
    drag = 0.5 * airDensity * area * C * speed**2
    lift = 0.5 * airDensity * area * CL * speed**2

    # Directions
    dragX = -drag * (velocityX / safeV)
    dragZ = -drag * (velocityZ / safeV)
    liftX = -lift * (velocityZ / safeV)
    liftZ = lift * (velocityX / safeV)

    # Update velocities
    velocityX += (dragX + liftX) * dt
    velocityZ += (-g + dragZ + liftZ) * dt

    # Update positions
    x += velocityX * dt
    z += velocityZ * dt
    t += dt

    trajectory.append((x, z))

    if x > 200:  # safety stop
        break

# ========== DRAW TRAJECTORY ==========
ctx.beginPath()
ctx.moveTo(50, 350)
scale_x = 800 / 300  # horizontal scaling
scale_z = 50         # vertical exaggeration

for px, pz in trajectory:
    canvas_x = 50 + px * scale_x
    canvas_y = 350 - pz * scale_z
    ctx.lineTo(canvas_x, canvas_y)

ctx.strokeStyle = "#0077ff"
ctx.lineWidth = 3
ctx.stroke()

# Ground line
ctx.beginPath()
ctx.moveTo(0, 350)
ctx.lineTo(800, 350)
ctx.strokeStyle = "#666"
ctx.lineWidth = 2
ctx.stroke()

# Ball marker
if trajectory:
    final_x, final_z = trajectory[-1]
    ctx.beginPath()
    ctx.arc(50 + final_x * scale_x, 350 - final_z * scale_z, 6, 0, 6.283)
    ctx.fillStyle = "#ff6600"
    ctx.fill()

# Output summary
dist_meters = trajectory[-1][0]
dist_yards = dist_meters * 1.0936
f"Total distance: {dist_yards:.1f} yards (t = {t:.2f}s)"
