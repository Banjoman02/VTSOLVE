"""__init__.py

Author: Monty Campbell
Created: A while ago
Description: Package Imports and Entry point."""

# Python Standard Imports
from math import degrees
import sys

# Third Party Imports
import numpy as np
from argparse import ArgumentParser

# Local Imports
from . import conversions, constants, pipeline, dynamics # NOTE: DO NOT DELETE OR YOU WILL BREAK THINGS
from .pipeline.data_importer import Measurement, loadData
from .constants import MU_SUN
from .pipeline.iod import iod
from .dynamics.propagation import intertial_cartesian_to_coe

def getParser() -> ArgumentParser:
    """Constructs the main CLI Parser.

    Returns:
        ArgumentParser: Argument parser.
    """
    parser = ArgumentParser(prog='Virginia Tech Space Object Location and Velocity Estimator',
                            usage='TODO: Add usage here.',
                            )
    
    parser.add_argument(
        'data_file',
        type=str,
        metavar='DATA',
        help='JSON file containing measurement data',
    )

    return parser


def runVTSOLVE(init_message_path:str) -> None:
    """Entry point.

    Args:
        init_message_path (str): Path to JSON file containing data.
    """
    # Load all the data
    print("="*80)
    print("Loading Measurements...")
    measurements:list[Measurement] = loadData(init_message_path)
    print(f"Loaded measurements!!!")
    for i,measurement in enumerate(measurements):
        print(f"Measurement {i}: {measurement}")

    # Put together Observer position matrix. NOTE: These are transposed per how Michael wrote the code.
    print("="*80)
    rog_array:np.ndarray = (np.array([measurement.observer_vector for measurement in measurements])).transpose() # Observer positions
    print("Calculated ROG Array!")
    print(f"{rog_array}")
    obs_array:np.ndarray = (np.array([measurement.rho_hat.vector for measurement in measurements])).transpose() # Unit vectors.
    print("Calculated boresight array!")
    print(f"{obs_array}")

    # Calculate time offsets.
    t1:float = (measurements[1].epoch - measurements[0].epoch).total_seconds()
    t3:float = (measurements[2].epoch - measurements[1].epoch).total_seconds()
    print("Calculated time offsets!")
    print(f"t1 = {t1} seconds")
    print(f"t3 = {t3} seconds")

    # Call IOD function.
    print("BEGINNING THE MAGIC!!!!")
    pos,vel = iod(obs_array,
                  rog_array,
                  t1,
                  t3,
                  mu=MU_SUN,
                  )
    print(f"Position and Velocity Calculated!\npos = {pos}\nvel={vel}")

    # State Vector to COE Conversions
    coe:tuple[float] = intertial_cartesian_to_coe(
        pos[0],
        pos[1],
        pos[2],
        vel[0],
        vel[1],
        vel[2],
        mu=MU_SUN,
    )
    print("="*80)
    print(f"Calculated Classical Orbital Elements!")
    print(f"SMA = {coe[0]} meters           ECC = {coe[1]}                INC = {degrees(coe[2])} deg")
    print(f"RAAN = {degrees(coe[3])} deg    ARGP = {degrees(coe[4])} deg  TANOM = {degrees(coe[5])} deg")

def main() -> None:
    parser = getParser()
    args = parser.parse_args()
    data_file:str = args.data_file

    runVTSOLVE(data_file)
    sys.exit(0)
