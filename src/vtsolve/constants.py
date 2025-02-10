"""constants.py

Author: Monty Campbell
Created: 2/3/2024
Description: Contains a few important constants lol."""

ECLIPTIC_INCLINATION: float = 0.4090926 # TODO: Cite this
"""``float``: Inclination of the ecliptic wrt the Earth's equator. Measured in radians."""

G:float = 6.67430e-11 # TODO: Cite this
"""``float``: Gravity constant in N * m^2 / kg^2."""

M_EARTH:float = 5.9722e24 # TODO: Cite this
"""``float``: Mass of the Earth in kilograms."""

M_SUN: float = 1.988416e30 # TODO: Cite this
"""``float``: Mass of the Sun in kilograms."""

MU_EARTH: float = G * M_EARTH
"""``float``: Gravity param for earth orbit."""

MU_SUN: float = G * M_SUN
"""``float``: Gravity param for the sun."""