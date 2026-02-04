import numpy as np
from datetime import datetime, timedelta
# Assumed date: February 3 2026

muSun = 1.32715*10**11  # km^3/s^2
semiMayAx = 2.68625*10**9   # km

halleyPeriod = 2*np.pi/np.sqrt(muSun)*semiMayAx**1.5

periodYears = halleyPeriod/(60*60*24*365.25)
periodLeftDays = (periodYears - int(periodYears))*365.25

print("Halley's comet period:", periodYears, "years")
print("Halley's comet period:", int(periodYears), "years and",
      periodLeftDays, "days")

timeInterv = (2026 - 1986)*365.25*24*60*60  # seconds

lastPerihelionDate = datetime(1986, 2, 9)
nextPerihelionDate = lastPerihelionDate + timedelta(seconds=halleyPeriod)
print("Next perihelion date (approx):", nextPerihelionDate.date())

halleyMeanMotion = np.sqrt(muSun/semiMayAx**3)  # rad
halleyMeanAnomaly = halleyMeanMotion*timeInterv     # rad


# Newton-Raphson solver for Kepler's equation: E - e sin(E) = M
def solve_kepler_eccentric_anomaly(M, e, tol=1e-10, max_iter=50):
    # Wrap M to [0, 2*pi) to aid convergence
    M = np.mod(M, 2*np.pi)
    # Initial guess
    E = M if e < 0.8 else np.pi
    for _ in range(max_iter):
        f = E - e*np.sin(E) - M
        fp = 1 - e*np.cos(E)
        dE = -f / fp
        E += dE
        if abs(dE) < tol:
            return E
    return E


eccentricity = 0.967
halleyEccentricAnomaly = solve_kepler_eccentric_anomaly(
    halleyMeanAnomaly, eccentricity
)

print("Halley's mean anomaly (rad):", halleyMeanAnomaly)
print("Halley's eccentric anomaly (rad):", halleyEccentricAnomaly)

# Find true anomaly and position for date
halleyTrueAnomaly = 2 * np.arctan2(
    np.sqrt(1 + eccentricity) * np.sin(halleyEccentricAnomaly / 2),
    np.sqrt(1 - eccentricity) * np.cos(halleyEccentricAnomaly / 2),
)
halleyRadius = semiMayAx * (1 - eccentricity * np.cos(halleyEccentricAnomaly))

halleyX = halleyRadius * np.cos(halleyTrueAnomaly)
halleyY = halleyRadius * np.sin(halleyTrueAnomaly)

print("Halley's true anomaly (rad):", halleyTrueAnomaly)
print("Halley's radius (A.U.):", halleyRadius*6.68459*10**(-9))
print("Halley's radius (km):", halleyRadius)
print("Halley's position in orbital plane (km):", (halleyX, halleyY))
