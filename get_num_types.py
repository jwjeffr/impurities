#!/usr/bin/env python

"""
small script for getting number of atom types from a LAMMPS data file
"""

import sys


def main():
    """
    get number of atom types from mass lines
    """

    num_types = 0
    with open(sys.argv[1], "r", encoding="utf8") as file:
        for line in file:
            if "mass" not in line:
                continue

            num_types += 1

    if not num_types:
        raise ValueError("num types not found")

    print(num_types)


if __name__ == "__main__":
    main()
