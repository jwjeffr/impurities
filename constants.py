"""
Module with some convenient constants
"""


class AtomColors:

    """
    Color constants for atoms
    """

    COBALT = (1.0, 0.4, 0.4)
    NICKEL = (0.4, 0.4, 1.0)
    CHROMIUM = (1.0, 1.0, 0.0)
    IRON = (1.0, 0.4, 1.0)
    MANGANESE = (0.4, 1.0, 0.2)
    ALUMINUM = (0.0, 1.0, 1.0)


class SystemColors:

    """
    Color constants for systems
    """

    CANTOR = "lightcoral"
    FEAL = "dodgerblue"


class SystemLineStyles:

    """
    Line style constants for systems
    """

    CANTOR = "-"
    FEAL = "--"


class Alignments:

    """
    Alignment flags for ovito overlays
    """

    TOP_LEFT = 33
    BOTTOM_LEFT = 65


SYSTEMS = ["cantor", "FeAl"]

NUM_FRAMES = 50

FINAL_STEP = 1_000_000

TYPE_MAPS = {
    "cantor": {1: "Co", 2: "Ni", 3: "Cr", 4: "Fe", 5: "Mn"},
    "FeAl": {1: "Fe", 2: "Al"},
}

NEAREST_NEIGHBORS = {"FCC": 12, "BCC": 8}
