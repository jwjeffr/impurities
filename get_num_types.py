#!/usr/bin/env python

import sys


def main():
    num_types = 0
    with open(sys.argv[1], "r") as file:
        for line in file:
            if "mass" not in line:
                continue

            num_types += 1

    if not num_types:
        raise ValueError("num types not found")

    print(num_types)


if __name__ == "__main__":
    main()
