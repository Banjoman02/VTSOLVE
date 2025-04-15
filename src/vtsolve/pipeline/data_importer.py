"""data_importer.py

Author: Monty Campbell
Created: 4/9/2025
Description: Imports necessary input data from a json file."""

# Python Standard Imports
from json import loads
from dataclasses import dataclass
from datetime import datetime
from math import radians

# Third Party Imports

# Local Imports
from ..conversions.coordinates import *
from ..dynamics.observer import Observer
from ..dynamics.propagation import OrbitalBodyPosition


@dataclass
class Measurement:
    """Encapsulates all data in a meausurement."""

    timestamp_iso:str
    """``str``: ISO Timestamp"""

    ra:float
    """``float``: Right ascension, in degrees."""

    dec:float
    """``float``: Declination, in degrees."""

    lat:float
    """``float``: Lattitude of the observer, in degrees."""

    lon:float
    """``float``: Longitude of the observer, in degrees."""

    @property
    def epoch(self) -> datetime:
        """``datetime``: Measurement epoch represented as a datetime object."""
        return datetime.fromisoformat(self.timestamp_iso)
    
    @property
    def celestial_coordinates(self) -> CelestialCoordinates:
        """``CelestialCoordinates``: Represents the meausurement as celestial coordinates."""
        return CelestialCoordinates(ra=radians(self.ra), dec=radians(self.dec))
    
    @property
    def rho_hat(self) -> EclipticEarthCenteredCoordinates:
        """``EclipticEarthCenteredCoordinates``: Expresses the meausurement as a line-of-sight vector."""
        return self.celestial_coordinates.toEarthCenteredEclipticUnit()
    
    @property
    def observer_vector(self) -> np.ndarray:
        """``np.ndarray``: Position from the sun pointing towards the observer's position on earth."""
        earth_orbit = OrbitalBodyPosition.earth()
        observer = Observer(radians(self.lat),
                            radians(self.lon),
                            )
        
        # Calculate earth's orbital position around the sun
        nu:float = earth_orbit.calcNu(self.epoch)
        r:float = earth_orbit.calcRad(nu)
        earth_pos:np.ndarray = earth_orbit.calcInertialCoords(nu, r)

        # Calculate observer's position wrt the earth in the solar frame.
        obs_pos = observer.toSolarInertial(self.epoch)

        # Calculate combined position and return
        return earth_pos + obs_pos

def loadData(path:str) -> list[Measurement]:
    """Loads data from input json.

    Args:
        path (str): Path to JSON file containing the data.

    Returns:
        list[Measurement]: List of all measurement objects.
    """
    with open(path, 'r') as fstream:
        data:list[dict] = loads(fstream)
    lat = float(data['lat'])
    lon = float(data['lon'])
    measurement_json:list[dict] = [measurement for measurement in data['measurements']]
    for measurement in measurement_json: # Populate JSON with observer position data.
        measurement['lat'] = lat
        measurement['lon'] = lon
    return [Measurement(**measurement) for measurement in measurement_json]
