"""coordinates.py

Author: Monty Campbell
Created: 2/3/24
Description: Contains some coordinate conversion functions lol.
"""

# Python Standard Imports
from math import sin, cos, pi

# Third Party Libraries
import numpy as np

# Local Imports
from ..constants import ECLIPTIC_INCLINATION

class SolarCenteredCoordinates:
    """Represents solar coordinates in the eclipltic frame, with x-hat pointing to J2000"""

class EclipticEarthCenteredCoordinates:
    """Represents cartesian coordinates in the ecliptic frame, but relative to the Earth's COM."""

    def __init__(self,
                 x: float,
                 y: float,
                 z: float,
                 ):
        """Constructs an instance.

        Args:
            x (float): X-coordinate, wrt the Earth, but in the ecliptic frame.
            y (float): Y-coordinate, wrt the Earth, but in the ecliptic frame.
            x (float): Z-coordinate, wrt the Earth, but in the ecliptic frame.
        """
        self._x:float = x
        """``float``: The x-coordinate, in the earth-centered ecliptic frame."""

        self._y:float = y
        """``float``: The y-coordinate, in the earth-centered ecliptic frame."""

        self._z:float = z
        """``float``: The z-coordinate, in teh earth-centered ecliptic frame."""

    @classmethod
    def fromVector(cls, vector:np.ndarray):
        """Constructs an instance from a position vector.

        Args:
            vector (np.ndarray): Position vector.
        """
        assert len(vector.tolist()) == 3,f"Position vector must be 3-dimensional. {len(vector.tolist())} dimensions given!"
        return cls(
            x=vector[0],
            y=vector[1],
            z=vector[2],
        )

    @property
    def vector(self) -> np.ndarray:
        """``np.ndarray``: Represents the coordinates as a cartesian vector."""
        return np.array([self._x, self._y, self._z])

class CelestialCoordinates:
    """Encapsulates the Celestial Coordinates"""

    def __init__(self,
                 ra:float,
                 dec:float,
                 ):
        """Constructs an instance.

        Args:
            ra (float): Right-ascension, in radians.
            dec (float): Declination, in radians.
        """
        self._ra:float = ra
        """``float``: Right-ascension, in radians."""

        self._dec:float = dec
        """``float``: Declination, in radians."""

    @property
    def vector(self) -> np.ndarray:
        """``np.ndarray``: Represents the set of coordinates as a cartesian vector. Note that the magnitude will be 1."""
        theta: float = self._ra
        phi:float = pi / 2 - self._dec

        x = cos(theta) * sin(phi)
        y = sin(theta) * sin(phi)
        z = cos(phi)

        return np.array([x, y, z])
    
    def toEarthCenteredEclipticUnit(self) -> EclipticEarthCenteredCoordinates:
        """Converts the celestial coordinates to earth centered ecliptic, but as a UNIT VECTOR!

        Returns:
            EclipticEarthCenteredCoordinates: Coordinates, represented as Earth Centered Ecliptic.
        """
        rot_matrix = np.array([1, 0, 0],
                              [0, cos(ECLIPTIC_INCLINATION), -1 * sin(ECLIPTIC_INCLINATION)],
                              [0, sin(ECLIPTIC_INCLINATION), cos(ECLIPTIC_INCLINATION)],
                              )
        converted = np.matmul(rot_matrix, self.vector)
        return EclipticEarthCenteredCoordinates.fromVector(converted)