import ovito
import numpy as np


def nearest_neighbor_topology_modifier(num_nearest_neighbors: int) -> callable:
    """
    Modifier for creating a topology from the N nearest neighbors
    """

    def wrapper(frame: int, data: ovito.data.DataCollection) -> None:
        finder = ovito.data.NearestNeighborFinder(num_nearest_neighbors, data)
        all_bonds, _ = finder.find_all()
        adjacency = np.zeros((max(all_bonds.shape), max(all_bonds.shape)))

        for i, row in enumerate(all_bonds):
            adjacency[i, row] = 1.0

        adjacency = np.triu(adjacency)
        bonds_container = data.particles_.create_bonds()
        bonds_container.create_property(
            "Topology", data=np.vstack(np.nonzero(adjacency)).T
        )

    return wrapper
