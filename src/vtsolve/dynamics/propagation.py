"""propagation.py
Author: Monty Campbell
Description: Does all the orbit propagation.
Created: 2/20/2025
"""

# Standard Library Imports
from math import *
import numpy as np
from datetime import datetime

# Local Imports
from ..constants import *
from .kepler import solveKepler

class OrbitalBodyPosition:
    
    def __init__(self,
                 sma: float,
                 ecc: float,
                 inc: float,
                 argp: float,
                 raan: float,
                 t0: datetime,
                 mean_lon: float,
                 mu:float = MU_SUN,
                 ):
        """Calculates the position of an orbital body.

        Args:
            sma (float): Semi major axis.
            ecc (float): Eccentricity.
            inc (float): Inclination.
            argp (float): Argument of Periapsis.
            raan (float): Right ascension of ascending node.
            t0 (float): Intial epoch.
            mean_lon (float): Mean Longitude.
            mu (float, Optional): Gravitational Parameter. Defaults to MU_SUN.
        """
        self.sma:float = sma
        self.ecc:float = ecc
        self.inc:float = inc
        self.argp:float = argp
        self.raan:float = raan
        self.t0:datetime = t0
        self.mean_lon:float = mean_lon
        
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
            J2000,
            EARTH_MEAN_LONG,
            mu=MU_SUN,
        )
        
    def calcNu(self, t:datetime) -> float:
        """Calculates true anomaly at a spot in time.

        Args:
            t (datetime): Future datetime.

        Returns:
            float: True anomaly at datetime.
        """
        delta_t:float = (t - self.t0).total_seconds() # Total time difference between inital time and input time
        mean_anom = self.mean_lon + sqrt(self.mu / (self.sma ** 3)) * delta_t - self.raan - self.argp # Mean anomaly
        ecc_anom = solveKepler(mean_anom, self.ecc) # Eccentric anomaly
        return 2 * atan2(sqrt((1 + self.ecc) / (1 - self.ecc)) * tan(ecc_anom / 2)) # Calculates and returns nu
    
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
                            ) -> np.ndarray:
        """Calculates coordinates in inertial frame.

        Args:
            nu (float): True anomaly.
            r (float): Orbital radius/

        Returns:
            np.ndarray: 3D Position vector containing coordinates in the inertial frame.
        """
        perifocal:np.ndarray = np.array([r * cos(nu), r * sin(nu), 0])
        R_3:np.ndarray = np.array([[cos(self.argp), sin(- self.argp), 0],
                                   [sin(self.argp), cos(self.argp), 0],
                                   [0, 0, 1],
                                   ])
        x_a:np.ndarray = np.matmul(R_3, perifocal)
        R_1:np.ndarray = np.array([
            [1, 0, 0],
            [0, cos(-self.inc), sin(-self.inc)],
            [0, -sin(-self.inc), cos(-self.inc)],
            ])
        x_b:np.ndarray = np.matmul(R_1, x_a)
        R_2:np.ndarray = np.array([
            [cos(-self.raan), 0, -sin(-self.raan)],
            [0, 1, 0],
            [sin(-self.raan), 0, cos(-self.raan)],
            ])
        x_i = np.matmul(R_2, x_b)
        return x_i # This is position in the inertial frame.
