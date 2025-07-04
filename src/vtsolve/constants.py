"""constants.py

Author: Monty Campbell
Created: 2/3/2024
Description: Contains a few important constants lol."""

# Python Standard Imports
from math import radians
from datetime import datetime, timezone

# ================ MISC ================

ECLIPTIC_INCLINATION: float = 0.4090926 # TODO: Cite this
"""``float``: Inclination of the ecliptic wrt the Earth's equator. Measured in radians."""

J2000: datetime = datetime(2000, 1, 1, 0, 0, 0).replace(tzinfo=timezone.utc)
"""``float``: Formal epoch definition of J2000."""

G:float = 6.67430e-11 # TODO: Cite this
"""``float``: Gravity constant in N * m^2 / kg^2."""

M_EARTH:float = 5.9722e24 # TODO: Cite this
"""``float``: Mass of the Earth in kilograms."""

M_SUN: float = 1.988416e30 # TODO: Cite this
"""``float``: Mass of the Sun in kilograms."""

MU_EARTH: float = G * M_EARTH
"""``float``: Gravity param for earth orbit."""

R_EARTH:float = 6378.137e3
"""``float``: Earth equatorial radius in meters."""

OMEGA_EARTH:float = 7.2921159e-5
"""``float``: Earth sidereal rotation angular velocity in rad/sec."""

EQUINOX_DT:datetime = datetime(2024, 3, 20, 3, 6, 20).replace(tzinfo=timezone.utc)
"""``datetime``: Datetime object of the last vernal equinox."""

# ================ EARTH ORBITAL PARAMS ================

MU_SUN: float = 1.327e20
"""``float``: Gravity param for the sun."""

EARTH_SMA: float = 149.598e9
"""``float``: Earth orbital sma around sun."""

EARTH_ECC:float = 0.01671022
"""``float``: Earth's orbital eccentricity."""

EARTH_INC: float = radians(0.00005)
"""``float``: Earth inclination in radians."""

EARTH_RAAN: float = radians(-11.26064)
"""``float``: Earth right ascension of ascending node."""

EARTH_ARGP: float = radians(102.94719)
"""``float``: Earth argument of periapsis."""

EARTH_M0: float = radians(357.51716)
"""``float``: Earth mean anomaly at epoch."""