"""data_importer.py

Author: Monty Campbell
Created: 4/9/2025
Description: Imports necessary input data from a json file."""

# Python Standard Imports
from json import loads
from dataclasses import dataclass
from datetime import datetime


@dataclass
class Measurement:
    """Encapsulates all data in a meausurement."""

    timestamp_iso:str
    """``str``: ISO Timestamp"""

    ra:float
    """``float``: Right ascension, in degrees."""

    dec:float
    """``float``: Declination, in degrees"""

    @property
    def epoch(self) -> datetime:
        """``datetime``: Measurement epoch represented as a datetime object."""
        return datetime.fromisoformat(self.timestamp_iso)
    

def loadData(path:str) -> list[Measurement]:
    """Loads data from input json.

    Args:
        path (str): Path to JSON file containing the data.

    Returns:
        list[Measurement]: List of all measurement objects.
    """
    with open(path, 'r') as fstream:
        data:list[dict] = loads(fstream)
    return [Measurement(**json_obj) for json_obj in data]
