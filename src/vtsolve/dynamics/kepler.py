"""kepler.py

Author: Monty Campbell
Created: 2/9/25
Description: Solves Kepler's equation.
"""

# Python Standard Imports
from math import * # LOL WILDCARD IMPORT HAHAHAHAH FUCK YOU PYTHON

ATOL:float = 1e-7
"""``float``: Absolute tolerance."""
MAX_ITER:int = 10000
"""``int``: Maximum iterations."""

class UnboundOrbitError(Exception):
    """Represents an exception for an unbound orbit."""
    pass

def _verifyBoundEcc(ecc: float) -> None:
    """Verifies that the eccentricity is for a bound orbit.

    Args:
        ecc (float): Eccentricity
    """
    if ecc < 0 or ecc >= 1:
        raise UnboundOrbitError(f"Orbit is not bound: ecc = {ecc}")

def _estErr(guess: float,
            ecc: float,
            ) -> float:
    """Estimates the error using Newton Raphson method.

    Args:
        guess (float): Guess for the eccentric anomaly.
        ecc (float): Eccentricity estimate.

    Returns:
        float: Estimated error.
    """
    # Calculate M as a function of E
    mean_anom: float = guess - ecc * sin(guess)
    # Calculate DM/DE
    d_mean_anom: float = 1 - ecc * cos(guess)
    return -1 * mean_anom / d_mean_anom

def solveKepler(init_guess: float,
                ecc: float,
                ) -> float:
    """Solves Kepler's Equation numerically.

    Args:
        init_guess (float): Initial guess for E. Typically M.
        ecc (float): Known eccentricity.

    Returns:
        float: Numerically approximated solution for E.
    """
    _verifyBoundEcc(ecc) # TODO: Verify if this actually fucking needs to be here.
    guess = init_guess + 0
    num_iter:int = 0
    while abs(init_guess - guess + ecc * sin(guess)) >= ATOL:
        num_iter += 1
        err_n: float = _estErr(guess, ecc)
        guess += err_n
        if num_iter >= MAX_ITER:
            break
    return guess
