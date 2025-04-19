"""propagation.py
Author: Monty Campbell
Description: Does all the orbit propagation.
Created: 2/20/2025
"""

# Standard Library Imports
from math import *
import numpy as np
from datetime import datetime, timezone

# Local Imports
from ..constants import *
from ..conversions.coordinates import R1, R2, R3
from .kepler import solveKepler

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
    raan:float = atan2(n_hat[1], n_hat[0])
    # argp
    argp:float = acos(np.dot(ecc_vec, n_hat) / e)
    if np.dot(ecc_vec, k_hat) < 0: # Quadrant check for z-component of the eccentricity vector.
        if debug:
            print("Quadrant check handled for argp.")
        argp = 2 * pi - argp
    # True Anomaly
    true_anom:float = acos(np.dot(ecc_vec, r) / (e * np.linalg.norm(r)))
    if np.dot(r, v) <= 0:
        true_anom = 2 * pi - true_anom # Quadrant check
    return (a, e, i, raan, argp, true_anom)

def main() -> None:
    """User interface for the calculations. Should not be used in the rest of the codebase."""
    r:list[float] = eval(input("Enter the position as a list with units in DU: "))
    v:list[float] = eval(input("Enter the velocity vector as a list with units of DU/TU: "))
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

class OrbitalBodyPosition:
    
    def __init__(self,
                 sma: float,
                 ecc: float,
                 inc: float,
                 argp: float,
                 raan: float,
                 M0:float,
                 mu:float = MU_SUN,
                 ):
        """Calculates the position of an orbital body.

        Args:
            sma (float): Semi major axis.
            ecc (float): Eccentricity.
            inc (float): Inclination.
            argp (float): Argument of Periapsis.
            raan (float): Right ascension of ascending node.
            M0 (float): Mean anomaly at epoch J2000.
            mu (float, Optional): Gravitational Parameter. Defaults to MU_SUN.
        """
        self.sma:float = sma
        self.ecc:float = ecc
        self.inc:float = inc
        self.argp:float = argp
        self.raan:float = raan
        self.m0:float = M0
        
        self.mu:float = mu

    @classmethod
    def earth(cls) -> "OrbitalBodyPosition":
        """Fast method for creating an instance that has all of the Earth's orbital elements pre-loaded.

        Returns:
            OrbitalBodyPosition: Constructed Earth Orbit.
        """
        return cls(
            EARTH_SMA,
            EARTH_ECC,
            EARTH_INC,
            EARTH_ARGP,
            EARTH_RAAN,
            EARTH_M0,
            mu=MU_SUN,
        )
        
    def calcNu(self, t:datetime) -> float:
        """Calculates true anomaly at a spot in time.

        Args:
            t (datetime): Future datetime.

        Returns:
            float: True anomaly at datetime.
        """
        delta_t:float = (t - J2000).total_seconds() # Total time difference between inital time and input time
        mean_anom = (self.m0 + sqrt(self.mu / (self.sma ** 3)) * delta_t) % (2 * pi)  # Mean anomaly
        ecc_anom = solveKepler(mean_anom, self.ecc) # Eccentric anomaly
        return 2 * atan2(sqrt((1 + self.ecc) / (1 - self.ecc)) * tan(ecc_anom / 2), 1) # Calculates and returns nu
    
    def calcRad(self, nu:float) -> float:
        """Calculates the orbital distance at a given true anomaly.

        Args:
            nu (float): True anomaly.

        Returns:
            float: Orbit radius.
        """
        return self.sma * (1 - self.ecc**2) / (1 + self.ecc * cos(nu))
    
    def calcInertialCoords(self,
                            nu:float,
                            r:float,
                            debug:bool = True,
                            ) -> np.ndarray:
        """Calculates coordinates in inertial frame.

        Args:
            nu (float): True anomaly.
            r (float): Orbital radius/

        Returns:
            np.ndarray: 3D Position vector containing coordinates in the inertial frame.
        """
        x_p:np.ndarray = np.array([r * cos(nu), r * sin(nu), 0]) # Perifocal coordinates
        x_a = np.matmul(R3(-self.argp), x_p) # Intermediate frame 1
        x_b = np.matmul(R1(-self.inc), x_a) # Intermediate frame 2
        x_i = np.matmul(R3(-self.raan), x_b) # Intertial Frame
        if debug:
            print(f"DEBUG x_i = {x_i}")
        return x_i # This is position in the inertial frame.
