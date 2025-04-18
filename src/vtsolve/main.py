#!/usr/bin/python

"""main.py

Author: Monty Campbell
Created: 4/17/25
Description: Main entry point and script for VT SOLVE"""

# Python Standard Imports
import sys

# Third Party Imports

# VT SOLVE
import vtsolve

def main() -> None:
    """Entry point"""
    parser = vtsolve.getParser()
    args = parser.parse_args()
    data_file:str = args.data_file

    vtsolve.main(data_file)
    sys.exit(0)
