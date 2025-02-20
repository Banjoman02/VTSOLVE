#!/usr/bin/pythonw
"""kepler.py

Author: Monty Campbell
Created: 2/7/2025
Description: Simple script for solving Kepler's equaution."""
# Standard Library Imports
from math import *
import sys

MAX_ITER:int = 10000
"""``int``: Maximum allowed number of iterations."""
ATOL:float = 1e-7
"""``float``: Tolerance."""

def calcError(Ecc_Anom: float,
              ecc: float,
              M: float,
              ) -> float:
    """Calculates the error for a given guess.

    Args:
        Ecc_Anom (float): Guess for eccentric anomaly
        ecc (float): Eccentricity.
        M (float): Mean Anomaly.

    Returns:
        float: Epsilon.
    """
    # Calculate f(x_n), or in this case M(E_n)
    mean_anom:float = Ecc_Anom - ecc * sin(Ecc_Anom) - M
    # Calculate f'(x_n), or in this case M'(E_n)
    d_mean_anom:float = 1 - ecc * cos(Ecc_Anom)
    # Calculate epsilon_n
    err = -1 * mean_anom / d_mean_anom
    return err

def solveKepler(init_guess:float,
                ecc: float,
                ) -> float:
    """Numerically solves Kepler's Equation.

    Args:
        init_guess (float): Initial guess for E. This is typically the Mean Anomaly, M.
        ecc (float): Orbit eccentricity.

    Returns:
        float: Best guess estimate for E.
    """
    E_n = init_guess
    for _i in range(MAX_ITER):
        E_last = E_n + 0
        E_n += calcError(E_n, ecc, init_guess)
        if abs(E_last - E_n) <= ATOL:
            print(f"Converged after {_i} iterations")
            return E_n
    raise RuntimeError(f"Failed to converge after {MAX_ITER} steps!")

def main() -> None:
    """What do you think this does lol."""
    args:list = sys.argv # When you write your own cli instead of using argparse...
    arg_format:str = "\nArgs must be as follows:\n python3 kepler.py [Mean Anomaly] [eccentricity]"
    if len(args) == 0: # The user gave us NOTHING
        # Tell them what a silly goose they are and make them give us things.
        print(f"No args given!\n{arg_format}")
        sys.exit(1)
    M = float(sys.argv[1])
    ecc = float(sys.argv[2])
    print("Solving for the Eccentric Anomaly.")
    try:
        E = solveKepler(M, ecc)
    except RecursionError:
        print(f"HAHAHAHAHAH LOL IT BROKE")
    print(f"Solved E = {E} rad = {degrees(E)} degrees!")
    sys.exit(0)

if __name__=='__main__':
    main()