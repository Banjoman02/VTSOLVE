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
from ressolve.image import Image, loadCommon
from argparse import ArgumentParser
import cv2
import numpy as np

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
        "-x",
        "--x-coordinate",
        dest="x",
        metavar="X",
        type=float,
        default=None,
        help="Asteroid Position X-Coordinate"
    )

    parser.add_argument(
        "-y",
        "--y-coordinate",
        dest="y",
        metavar="Y",
        type=float,
        default=None,
        help="Asteroid Position Y-Coordinate",
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

def drawCircle(image:Image,
               x:int,
               y:int,
               radius=50,
               color=255,
               thickness=2) -> None:
    """_summary_

    Args:
        image (Image): _description_
        x (int): _description_
        y (int): _description_
        radius (int, optional): _description_. Defaults to 50.
        color (int, optional): _description_. Defaults to 255.
        thickness (int, optional): _description_. Defaults to 2.
    """
    image_with_circle = image.img_arr.copy()
    center = (int(x), int(y))
    cv2.circle(image_with_circle, center, radius, color, thickness)
    image.img_arr = image_with_circle

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
    mirrored:bool = args.mirrored
    x:float | None = args.x
    y:float | None = args.y

    # Open the image
    _img_arr = loadCommon(img_path)
    image = Image.fromImageArray(_img_arr)
    if mirrored:
        image.flipHorizontal()
        image.saveToFITS(save_path=img_path)
    
    # Plate-solve the image
    solver = AstrometrySolver(image, img_path=img_path)
    solution:PlateSolution = solver.solveField()

    # Print some stuff for the solved RA/DEC
    print(f"Imaged Solved = {solution.solved}\n")
    print(f"Center RA = {solution.center_ra}. Center DEC = {solution.center_dec}")
    print(f"Image is rotated {solution.center_rotation} degrees E of N.")

    if x is not None and y is not None:
        coords = solver.getPixelPosition(x, y)
        print(f"Solved coords at (x, y) = {coords}")

if __name__=='__main__':
    main()