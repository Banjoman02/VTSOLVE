"""Put something here later"""

# Python Standard Imports
from math import *
import numpy as np
from datetime import datetime

# Local Imports
from ..constants import R_EARTH, OMEGA_EARTH, EQUINIOX_DT, ECLIPTIC_INCLINATION

class Observer:
    
    def __init__(self,
                 lat: float,
                 lon: float,
                 rad: float = R_EARTH,
                 ):
        """_summary_

        Args:
            lat (float): _description_
            lon (float): _description_
            rad (float, optional): _description_. Defaults to R_EARTH.
        """
        self.lat:float = lat
        self.lon:float = lon
        self.radius:float = rad
        
        self.gamma = pi / 2 - self.lat
    
    def toECI(self, t:datetime) -> np.ndarray:
        """_summary_

        Returns:
            np.ndarray: _description_
        """
        ecef:np.ndarray =  np.array([
            self.radius * cos(self.lon) * sin(self.gamma),
            self.radius * sin(self.lon) * sin(self.gamma),
            self.radius * cos(self.gamma),
        ])
        delta_t = (t - EQUINIOX_DT).total_seconds()
        theta = (OMEGA_EARTH * delta_t) % (2 * pi)
        R_3 = np.array([
            [cos(theta), sin(theta), 0],
            [-sin(theta), cos(theta), 0],
            [0, 0, 1],
            ])
        return np.matmul(R_3, ecef)
    
    def toSolarInertial(self,
                        t:datetime,
                        ) -> np.ndarray:
        """_summary_

        Args:
            t (datetime): _description_

        Returns:
            np.ndarray: _description_
        """
        eci:np.ndarray = self.toECI(t)
        rot_matrix = np.array([1, 0, 0],
                              [0, cos(ECLIPTIC_INCLINATION), -1 * sin(ECLIPTIC_INCLINATION)],
                              [0, sin(ECLIPTIC_INCLINATION), cos(ECLIPTIC_INCLINATION)],
                              )
        return np.matmul(rot_matrix, eci)    
        