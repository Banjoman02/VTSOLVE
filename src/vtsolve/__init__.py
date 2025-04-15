"""__init__.py

Author: Monty Campbell
Created: A while ago
Description: Package Imports and Entry point."""

# Python Standard Imports
from datetime import datetime, timedelta

# Third Party Imports
import numpy as np

# Local Imports
from . import conversions, constants, pipeline, dynamics # NOTE: DO NOT DELETE OR YOU WILL BREAK THINGS
from .pipeline.data_importer import Measurement, loadData
from .constants import MU_SUN
from .pipeline.iod import iod

def main(init_message_path:str) -> None:
    """Entry point.

    Args:
        init_message_path (str): Path to JSON file containing data.
    """
    # Load all the data
    print("Loading Measurements...")
    measurements:list[Measurement] = loadData(init_message_path)
    print("Loaded measurements!!!")

    # Put together Observer position matrix. NOTE: These are transposed per how Michael wrote the code.
    rog_array:np.ndarray = (np.array([measurement.observer_vector for measurement in measurements])).transpose() # Observer positions
    print("Calculated ROG Array!")
    obs_array:np.ndarray = (np.array([measurement.rho_hat for measurement in measurements])).transpose() # Unit vectors.
    print("Calculated boresight array!")

    # Calculate time offsets.
    t1:float = (measurements[1].epoch - measurements[0].epoch).total_seconds()
    t3:float = (measurements[2].epoch - measurements[1].epoch).total_seconds()
    print("Calculated time offsets!")

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
