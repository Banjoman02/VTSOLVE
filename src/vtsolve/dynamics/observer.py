"""Put something here later"""

# Python Standard Imports
from math import *
import numpy as np
from datetime import datetime, timezone

# Local Imports
from ..constants import R_EARTH, OMEGA_EARTH, EQUINOX_DT, ECLIPTIC_INCLINATION

class Observer:
    
    def __init__(self,
                 lat: float,
                 lon: float,
                 rad: float = R_EARTH,
                 ):
        """Constructs an instance.

        Args:
            lat (float): Geographic latitude of the observer.
            lon (float): Geographic longitude of the observer.
            rad (float, optional): Equatorial radius of the Earth. Defaults to R_EARTH.
        """
        self.lat:float = lat
        self.lon:float = lon
        self.radius:float = rad
        
        self.gamma = pi / 2 - self.lat
    
    def toECI(self, t:datetime) -> np.ndarray:
        """Computes the relative position of the observer in the ECI frame.

        Returns:
            np.ndarray: 3D Position vector in the ECI frame of the observer.
        """
        ecef:np.ndarray =  np.array([
            self.radius * cos(self.lon) * sin(self.gamma),
            self.radius * sin(self.lon) * sin(self.gamma),
            self.radius * cos(self.gamma),
        ])
        t.replace(tzinfo=timezone.utc)
        delta_t = (t - EQUINOX_DT).total_seconds()
        theta = (OMEGA_EARTH * delta_t) % (2 * pi)
        R_3 = np.array([
            [cos(theta), sin(theta), 0],
            [-sin(theta), cos(theta), 0],
            [0, 0, 1],
            ])
        return np.matmul(R_3, ecef)
    
    def toSolarInertial(self,
                        t:datetime,
                        debug:bool = False, 
                        ) -> np.ndarray:
        """Computes the coordinates in the solar inertial frame.

        Args:
            t (datetime): Epoch.

        Returns:
            np.ndarray: 3D Relative position vector in the solar frame.
        """
        eci:np.ndarray = self.toECI(t)
        rot_matrix = np.array([[1, 0, 0],
                              [0, cos(ECLIPTIC_INCLINATION), -1 * sin(ECLIPTIC_INCLINATION)],
                              [0, sin(ECLIPTIC_INCLINATION), cos(ECLIPTIC_INCLINATION)],],
                              )
        sci = np.matmul(rot_matrix, eci)
        if debug:
            print(f"Observer POS Debug: eci = {eci}, sci = {sci} ALL RELATIVE TO EARTH")
        return sci
        