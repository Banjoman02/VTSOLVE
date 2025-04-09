#!/usr/bin/bash

"""get_asteroid_pos.py
Author: Monty Campbell
Created: 4/8/2025
Description: Simple script that plate-solves an image and attempts to identify an object closest to a set of 
reported coordinates."""

# Python Standard Imports
import math

# Third Party Imports
from ressolve.platesolve.solution import PlateSolution
from ressolve.platesolve.astrometry_solver import AstrometrySolver
from ressolve.image import Image, loadFITS
from argparse import ArgumentParser

# Local Imports


def getCLI() -> ArgumentParser:
    """Builds the argument parser.

    Returns:
        ArgumentParser: Argument parser.
    """
    parser = ArgumentParser()

    parser.add_argument(
        "target_file",
        metavar="TARGET_FILE",
        type=str,
        help="Path to target file.",
    )

    parser.add_argument(
        "-ra",
        "--right-ascension",
        dest="ra",
        metavar="RA",
        type=float,
        default=None,
        help="Predicted right ascension, in degrees."
    )

    parser.add_argument(
        "-dec",
        "--declination",
        dest="dec",
        metavar="DEC",
        type=float,
        default=None,
        help="Predicted declination, in degrees",
    )

    parser.add_argument(
        "-m",
        "--mirrored",
        dest="mirrored",
        metavar="MIRRORED",
        default=False,
        type=bool,
        help="Optional flag if the image is mirrored. Will flip the image array before attempting to plate-solve."
    )

    return parser

def angular_distance(ra1_deg, dec1_deg, ra2_deg, dec2_deg):
    """
    Calculate the angular distance between two RA/Dec coordinates in degrees.
    
    Parameters:
    - ra1_deg, dec1_deg: coordinates of the first object (degrees)
    - ra2_deg, dec2_deg: coordinates of the second object (degrees)
    
    Returns:
    - Angular distance in degrees
    """
    # Convert degrees to radians
    ra1 = math.radians(ra1_deg)
    dec1 = math.radians(dec1_deg)
    ra2 = math.radians(ra2_deg)
    dec2 = math.radians(dec2_deg)
    
    # Spherical law of cosines
    cos_angle = math.sin(dec1) * math.sin(dec2) + math.cos(dec1) * math.cos(dec2) * math.cos(ra1 - ra2)
    # Clamp to valid range in case of numerical errors
    cos_angle = min(1.0, max(-1.0, cos_angle))
    angle_rad = math.acos(cos_angle)
    
    # Convert radians back to degrees
    angle_deg = math.degrees(angle_rad)
    return angle_deg

def main() -> None:
    """Entry Point"""
    parser:ArgumentParser = getCLI()
    args = parser.parse_args()

    print("Solving image!!!")

    # Pull required information from CLI
    img_path:str = args.target_file
    ra:float | None = args.ra
    dec:float | None = args.dec
    mirrored:bool = args.mirrored

    # Open the image
    _img_arr = loadFITS(img_path)
    image = Image.fromImageArray(_img_arr)
    if mirrored:
        image.flipHorizontal()
        image.saveToFITS(save_path=img_path)
    
    # Plate-solve the image
    solver = AstrometrySolver(image, img_path=img_path)
    solution:PlateSolution = solver.solveField()

    # Print some stuff for the solved RA/DEC
    print("Imaged Solved!\n")
    print(f"Center RA = {solution.center_ra}. Center DEC = {solution.center_dec}")
    print(f"Image is rotated {solution.center_rotation} degrees E of N.")

    if ra is not None and dec is not None:
        dist = angular_distance(ra, dec, solution.center_ra, solution.center_dec)
        print(f"Image Center is {dist} degrees away from center coordinates!")

if __name__=='__main__':
    main()