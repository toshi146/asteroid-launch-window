
"""
Asteroid Rendezvous Launch Window Tool
--------------------------------------
Compatible with poliastro 0.7.0
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

from astropy import units as u
from astropy.time import Time

from poliastro.bodies import Sun, Earth
from poliastro.twobody import Orbit
from poliastro.iod import lambert as lambert_iod  


# ------------------------- input helpers ----------------------------

def ask_float(prompt):
    while True:
        try:
            return float(input(prompt).strip())
        except ValueError:
            print("  ! Enter a number")


def ask_int(prompt, min_value=None):
    while True:
        try:
            v = int(input(prompt).strip())
            if min_value is not None and v < min_value:
                print(f"  ! Must be >= {min_value}")
                continue
            return v
        except ValueError:
            print("  ! Enter an integer")


def ask_time(prompt):
    while True:
        try:
            return Time(input(prompt).strip(), scale="tdb")
        except Exception:
            print("  ! Invalid date (e.g. 2030-01-01 00:00)")


# -------------------- anomaly conversion ----------------------------

def mean_to_true(M_deg, e):
    M = np.deg2rad(M_deg)
    e = float(e)

    E = M if e < 0.8 else np.pi
    for _ in range(50):
        dE = -(E - e*np.sin(E) - M) / (1 - e*np.cos(E))
        E += dE
        if abs(dE) < 1e-12:
            break

    sinv = np.sqrt(1 - e*e) * np.sin(E) / (1 - e*np.cos(E))
    cosv = (np.cos(E) - e) / (1 - e*np.cos(E))
    return np.arctan2(sinv, cosv)


# ---------------------------- logic ---------------------------------

def build_asteroid_orbit():
    print("\n--- Asteroid Elements ---")
    a = ask_float("a [AU]: ")
    e = ask_float("e [-]: ")
    i = ask_float("i [deg]: ")
    raan = ask_float("RAAN Ω [deg]: ")
    argp = ask_float("ω [deg]: ")
    M_deg = ask_float("M [deg]: ")
    epoch = ask_time("Elements epoch: ")

    nu = mean_to_true(M_deg, e) * u.rad

    return Orbit.from_classical(
        Sun, a*u.AU, e*u.one, i*u.deg, raan*u.deg, argp*u.deg, nu, epoch=epoch
    )


def get_params():
    print("\n--- Mission Parameters ---")
    start = ask_time("Launch window START: ")
    end = ask_time("Launch window END: ")

    if end <= start:
        print("End must be after start")
        sys.exit(1)

    tof_days = ask_float("Time of flight [days]: ")
    samples = ask_int("Number of departures to sample: ", 2)

    return start, end, tof_days, samples


def compute_curve(asteroid, dep_times, tof_days):
    n = len(dep_times)
    dv = np.full(n, np.nan)

    print("\n--- Running Lambert transfers ---\n")

    # mu as Quantity, not float
    mu = Sun.k.to(u.km**3 / u.s**2)

    for k, dep in enumerate(dep_times):
        arr = dep + tof_days * u.day

        # Earth at departure
        try:
            earth = Orbit.from_body_ephem(Earth, dep)
        except Exception as e:
            print(f"[{k+1}/{n}] Earth ephem failed: {e}")
            continue

        # Asteroid at arrival
        try:
            dt = (arr - asteroid.epoch).to(u.s)
            ast = asteroid.propagate(dt)
        except Exception as e:
            print(f"[{k+1}/{n}] Asteroid propagation failed: {e}")
            continue

        try:
            # Positions and velocities as Quantities
            r1 = earth.r.to(u.km)
            v1 = earth.v.to(u.km / u.s)
            r2 = ast.r.to(u.km)
            tof = tof_days * u.day

        
            sols = list(lambert_iod(mu, r1, r2, tof))
            if len(sols) == 0:
                raise RuntimeError("No Lambert solution returned")
            v_dep, v_arr = sols[0]

      
            dv1 = np.linalg.norm((v_dep - v1).to(u.km / u.s).value)
            dv2 = np.linalg.norm((v_arr - ast.v.to(u.km / u.s)).value)
            dv[k] = dv1 + dv2


            print(f"[{k+1:2d}/{n}] Δv = {dv[k]:.3f} km/s")

        except Exception as e:
            print(f"[{k+1}/{n}] Lambert failed: {e}")
            continue

    return dv


def plot_curve(times, dv):
    mask = ~np.isnan(dv)
    if not np.any(mask):
        print("\nNo valid Lambert data to plot.")
        return

    tx = [t.to_datetime() for t in times[mask]]

    plt.figure(figsize=(10, 5))
    plt.plot(tx, dv[mask], marker="o")
    plt.grid(True)
    plt.xlabel("Departure date")
    plt.ylabel("Total Δv [km/s]")
    plt.title("Launch Window Δv Curve")
    plt.tight_layout()
    plt.show()


def main():
    print("===========================================")
    print("  Asteroid Rendezvous Launch Window Tool")
    print("  (poliastro 0.7.0 compatible)")
    print("===========================================")

    asteroid = build_asteroid_orbit()
    start, end, tof, N = get_params()

    dep_jd = np.linspace(start.jd, end.jd, N)
    dep_times = Time(dep_jd, format="jd", scale="tdb")

    dv = compute_curve(asteroid, dep_times, tof)

    if np.all(np.isnan(dv)):
        print("\nNo valid Lambert solutions.")
    else:
        idx = np.nanargmin(dv)
        best = dep_times[idx]
        bestdv = dv[idx]
        arr = best + tof * u.day

        print("\n===== BEST SOLUTION =====")
        print("Departure:", best.iso)
        print("Arrival:  ", arr.iso)
        print(f"Min Δv:   {bestdv:.3f} km/s")

    plot_curve(dep_times, dv)
    print("\nDone.\n")


if __name__ == "__main__":
    main()
