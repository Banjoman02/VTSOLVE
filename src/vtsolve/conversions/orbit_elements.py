"""orbital_elements.py
Author: Monty Campbell
Created: 2/23/2025
Description: Does conversion between Classical Orbital Elements and Inertial coordinates."""

# Standard Library Imports
from math import *

# Third Party Imports
import numpy as np

# Local Imports

def intertial_cartesian_to_coe(
        x:float,
        y:float,
        z:float,
        vx:float,
        vy:float,
        vz:float,
        mu:float,
        debug:bool=False,
    ) -> tuple[float]:
    """Converts position and velocity into classical orbital elements.

    Args:
        x (float): X-coordinate, in meters.
        y (float): Y-coordinate, in meters.
        z (float): Z-coordinate, in meters.
        vx (float): X-velocity, in meters/second.
        vy (float): Y-velocity, in meters/second.
        vz (float): Z-velocity, in meters/second.
        mu (float, optional): Gravity Parameter, in units of m^3s^-2. Defaults to MU_SUN.

    Returns:
        tuple[float]: 
            (semi-major axis (a) in meters,
            eccentricity,
            inclination,
            raan,
            argp,
            true_anom,
            )
    """
    # Get position and velocity in vector form.
    r = np.array([x, y ,z])
    v = np.array([vx, vy, vz])
    # Unit vectors for the frame
    i_hat = np.array([1, 0, 0]) # i-hat vector
    j_hat = np.array([0, 1, 0]) # j-hat vector
    k_hat = np.array([0, 0, 1]) # k-hat vector
    # a, e, and i (the easy ones lol)
    epsilon:float = 0.5 * np.linalg.norm(v) ** 2 - mu / np.linalg.norm(r) # Specific Mechanical Energy
    h:np.ndarray = np.cross(r, v) # Specific angular momentum
    n:np.ndarray = np.cross(k_hat, h)
    n_hat:np.ndarray = n / np.linalg.norm(n) # Unit vector for the line of nodes.
    a:float = -1 * mu / (2 * epsilon) # Semi major axis
    ecc_vec:np.ndarray = 1 / mu * (np.cross(v, h) - (mu * r) / np.linalg.norm(r)) # Calculate eccentricity vector
    e:float = np.linalg.norm(ecc_vec) # Eccentricity as a scalar parameter
    i = acos(np.dot(k_hat, h) / np.linalg.norm(h)) # Inclination
    # RAAN
    raan:float = acos(np.dot(i_hat, h) / np.linalg.norm(h))
    if np.dot(j_hat, h) < 0: # Quadrant check for y-component of h
        raan = 2 * pi - raan
    # argp
    argp:float = acos(np.dot(ecc_vec, n_hat) / e)
    if np.dot(ecc_vec, k_hat) < 0: # Quadrant check for z-component of the eccentricity vector.
        argp = 2 * pi - argp
    # True Anomaly
    true_anom:float = acos(np.dot(ecc_vec, r) / (e * np.linalg.norm(r))) # No quadrant check needed here
    return (a, e, i, raan, argp, true_anom)

def main() -> None:
    """User interface for the calculations. Should not be used in the rest of the codebase."""
    r:list[float] = eval(input("Enter the position as a list with units in DU: "))
    v:list[float] = eval(input("Enter the velocity vector as a list with units of DU/TU:" ))
    mu:float = eval(input("Enter your gravity parameter in units if DU^3 / TU^2: "))
    coe = intertial_cartesian_to_coe(
        r[0],
        r[1],
        r[2],
        v[0],
        v[1],
        v[2],
        mu=mu,
    )
    print(f"a = {coe[0]}, e = {coe[1]}, i = {coe[2]}")
    print(f"raan = {coe[3]}, argp = {coe[4]}, true_anom = {coe[5]}")

if __name__=='__main__':
    main()