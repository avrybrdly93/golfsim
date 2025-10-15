import math

# Sample input values (you can update or hook to HTML later)
firstRun = True
spinRate = 3000  # RPM
phideg = 15
initialVelocity = 70  # m/s
airDensity = 1.225  # kg/m^3
airViscosity = 1.81e-5
radius = 0.02135  # m (e.g., golf ball)
dt = 0.01  # timestep
hasDimples = True

# Initial velocities (3D)
velocityY = 0
wx, wy, wz = 0, 0, 0  # Spin components
x, y, z = 0, 0, 1  # Starting position
t = 0
count = 0

# Physics calculations
if firstRun:
    omega = spinRate * 2 * math.pi / 60.0
    print("This is the first time!")
    phirad = math.pi * phideg / 180
    velocity = initialVelocity
    velocityX = initialVelocity * math.cos(phirad)
    velocityZ = initialVelocity * math.sin(phirad)
    print(f"First evolution, spinrate, launch angle, and initial velocity: {omega:.2f}, {velocity:.2f}, {phideg}")
    firstRun = False

# Constants
EPS_V = 0.1
S_MAX = 0.4
CL_MIN = 0.0
CL_MAX = 1.0
decayRate = 4.5

speed = math.sqrt(velocityX**2 + velocityY**2 + velocityZ**2)
safeV = max(speed, EPS_V)

reynoldsNumber = (airDensity * speed * radius * 2) / airViscosity

S = (omega * radius) / safeV
omega *= math.exp(-dt / decayRate)
print(f"Omega after decay: {omega:.2f} rad/s")

if hasDimples:
    CL = 1.99 * S - 3.25 * S**2
    if reynoldsNumber > 175000:
        C = min(1.91e-11 * reynoldsNumber**2 - 5.40e-6 * reynoldsNumber + 0.56, 0.21)
    elif reynoldsNumber > 80000:
        C = 1.91e-11 * reynoldsNumber**2 - 5.40e-6 * reynoldsNumber + 0.56
    else:
        print(f"Low-speed regime: velocity = {speed:.2f}")
        C = 0.5  # Fallback drag coefficient
else:
    C = 0.5
    CL = -0.38 * S**2 + 0.43 * S

area = math.pi * radius**2
liftMag = 0.5 * airDensity * area * abs(safeV) * safeV * CL

# Cross product omega x velocity
crossX = wy * velocityZ - wz * velocityY
crossY = wz * velocityX - wx * velocityZ
crossZ = wx * velocityY - wy * velocityX
crossMag = math.sqrt(crossX**2 + crossY**2 + crossZ**2)

signZ = -1 if velocityZ < 0 else 1

if crossMag < 1e-8:
    magnusForceX = magnusForceY = magnusForceZ = 0.0
else:
    magnusForceX = liftMag * (crossX / crossMag)
    magnusForceY = -liftMag * (crossY / crossMag)
    magnusForceZ = -liftMag * (crossZ / crossMag)

dragMag = 0.5 * airDensity * area * C * velocity**2
dragForceX = dragMag * (velocityX / velocity)
dragForceY = dragMag * (velocityY / velocity)
dragForceZ = dragMag * (velocityZ / velocity)

TotalDist = x * 1.0936
if z <= 0.05 and x > 5 and count == 0:
    TotalDist = round(100 * TotalDist) / 100
    PlaceX = x
    PlaceY = y
    TotalDistString = f"{TotalDist} yards"
    print(f"Total Distance Traveled: {round(TotalDist)}")
    print(f"Time: {t}")
    count += 1

# prints out coefficients and magnus force
print(f"Drag Coefficient (C): {C:.3f}")
print(f"Lift Coefficient (CL): {CL:.3f}")
print(f"Magnus Force X: {magnusForceX:.3f}")
