"""
Module for performing site statistics calculations, including chemical potential,
vacancy concentration, vacancy formation energy, and vacancy formation volume.
"""


from dataclasses import dataclass
from typing import Annotated, Tuple
from itertools import combinations
import numpy as np
from numpy.typing import ArrayLike


# some convenient variable names to define shapes in type annotations
NUM_TYPES = ...
NUM_SITES = ...
NUM_TEMPERATURES = ...


@dataclass
class FormationCalculator:

    """
    class for computing vacancy concentration and formation characteristics
    """

    energetics_data: Annotated[ArrayLike, (NUM_TYPES, NUM_SITES)]
    volumetrics_data: Annotated[ArrayLike, (NUM_TYPES, NUM_SITES)]

    def __post_init__(self):
        """
        checks that data arrays have the same shape, and vectorizes concentration calculation
        and formation quantity calculation at multiple temperatures
        """

        assert self.energetics_data.shape == self.volumetrics_data.shape

        # vectorize concentration and formation quantity methods
        # excluding all arguments except beta
        self.concentration_vectorized: callable = np.vectorize(
            self.get_concentration, excluded="self"
        )
        self.formation_quantity_vectorized: callable = np.vectorize(
            self.get_formation_quantity, excluded=["self", "quantity"]
        )

    def get_concentration(self, beta: float) -> float:
        """
        get vacancy concentration at a given temperature
        :param beta: inverse temperature
        :return: concentration at beta
        """

        # get boltzmann factors for each type-site pair
        boltzmann_factors = np.exp(beta * self.energetics_data)

        # sum over elements per site
        sum_of_boltzmann_factors = np.sum(boltzmann_factors, axis=0)

        # average 1 / (1 + sum)
        local_probabilities = 1.0 / (1.0 + sum_of_boltzmann_factors)

        return np.mean(local_probabilities)

    def get_formation_quantity(
        self, beta: float, quantity: Annotated[ArrayLike, (NUM_TYPES, NUM_SITES)]
    ) -> float:
        """
        general method for computing a formation quantity at a given temperature
        :param beta: inverse temperature
        :param quantity: thermodynamic formation quantity of interest
        :return: given formation quantity at beta
        """

        # get boltzmann factors for each type-site pair
        boltzmann_factors = np.exp(beta * self.energetics_data)

        # perform sum in numerator of equation defining global formation energy
        numerator = np.sum(quantity * boltzmann_factors, axis=0)

        # calculate denominator in the same equation
        sum_of_boltzmann_factors = np.sum(boltzmann_factors, axis=0)
        denominator = (1.0 + sum_of_boltzmann_factors) ** 2

        # perform thermodynamic average in the same equation
        thermodynamic_average = np.mean(numerator / denominator)

        # return average over concentration
        return thermodynamic_average / self.concentration_vectorized(beta)

    def get_formation_energy(self, beta: float) -> float:
        """
        method for getting formation energy at a given temperature
        :param beta: inverse temperature
        :return: formation energy at beta
        """

        return self.formation_quantity_vectorized(beta, quantity=self.energetics_data)

    def get_formation_volume(self, beta: float) -> float:
        """
        method for getting formation volume at a given temperature
        :param beta: inverse temperature
        :return: formation volume at beta
        """

        return self.formation_quantity_vectorized(beta, quantity=self.volumetrics_data)


def get_chemical_potentials(
    types: Annotated[ArrayLike, NUM_TYPES],
    occupying_energies: Annotated[ArrayLike, (NUM_TYPES, NUM_SITES)],
    concentrations: Annotated[ArrayLike, NUM_TYPES],
    enthalpy_per_atom: float,
) -> Annotated[ArrayLike, NUM_TYPES]:
    """
    method for getting chemical potentials using site statistics
    :param types: array of atom types, shape is (number of types,)
    :param occupying_energies: array of occupying energies
        shape is (number of types, number of sites)
    :param concentrations: array of concentrations in at. %, shape is (number of types,)
    :param enthalpy_per_atom: enthalpy per atom of equilibrated configuration
    :return: chemical potential of each type, shape is (number of types,)
    """

    assert types.shape == concentrations.shape == (occupying_energies.shape[0],)

    # initialize list of unique type pairs alpha != alpha'
    type_pairs = list(combinations(types, r=2))
    num_pairs = len(type_pairs)
    num_types = len(types)

    # initialize coefficient matrix to solve Ax = b
    # chemical potentials are unknown
    coefficient_matrix = np.zeros((num_pairs + 1, num_types))
    b = np.zeros(num_pairs + 1)

    # populate coefficient matrix
    for index, (first_type, second_type) in enumerate(type_pairs):
        coefficient_matrix[index, first_type] = 1.0
        coefficient_matrix[index, second_type] = -1.0
        b[index] = np.mean(
            occupying_energies[first_type, :] - occupying_energies[second_type, :]
        )

    # populate last members of coefficient matrix
    coefficient_matrix[num_pairs, :] = concentrations
    b[num_pairs] = enthalpy_per_atom

    # solve and return the least squares solution
    return np.linalg.lstsq(coefficient_matrix, b, rcond=None)[0]


def get_formation_array(
    occupying: Annotated[ArrayLike, (NUM_TYPES, NUM_SITES)],
    vacant: Annotated[ArrayLike, NUM_SITES],
    chemical_potentials: Annotated[ArrayLike, NUM_TYPES] = None,
) -> Annotated[ArrayLike, (NUM_TYPES, NUM_SITES)]:
    """
    helper method for getting formation array, i.e. penalties for each site and each type
    :param occupying: occupying array, shape (number of types, number of sites)
    :param vacant: vacant array, shape (number of sites,)
    :param chemical_potentials: optional chemical potential array, shape (number of types,)
    :return: formation array, shape (number of types, number of sites)
    """

    assert vacant.shape == (occupying.shape[1],)

    # get formation array by subtracting occupying values from vacant volumes
    num_types, num_sites = occupying.shape
    formation_array = np.vstack([vacant] * num_types) - occupying

    # if chemical potentials provided, add them to formation array
    # extra energetic penalty with chemical potentials
    if chemical_potentials is not None:
        assert chemical_potentials.shape == (occupying.shape[0],)
        formation_array += np.vstack([chemical_potentials] * num_sites).T

    return formation_array


def get_vacancy_characteristics(
    vacant_energies: Annotated[ArrayLike, NUM_SITES],
    occupying_energies: Annotated[ArrayLike, (NUM_TYPES, NUM_SITES)],
    vacant_volumes: Annotated[ArrayLike, NUM_SITES],
    occupying_volumes: Annotated[ArrayLike, (NUM_TYPES, NUM_SITES)],
    types: Annotated[ArrayLike, NUM_TYPES],
    concentrations: Annotated[ArrayLike, NUM_TYPES],
    enthalpy_per_atom: float,
    beta_vals: Annotated[ArrayLike, NUM_TEMPERATURES],
) -> Tuple[
    Annotated[ArrayLike, NUM_TEMPERATURES],
    Annotated[ArrayLike, NUM_TEMPERATURES],
    Annotated[ArrayLike, NUM_TEMPERATURES],
]:
    """
    helper method for getting vacancy characteristics
        (concentration and formation quantities) at multiple temperatures
    :param vacant_energies: array of vacant energies, shape is (number of sites,)
    :param occupying_energies: array of occupying energies
        shape is (number of types, number of sites)
    :param vacant_volumes: array of vacant volumes, shape is (number of sites,)
    :param occupying_volumes: array of occupying volumes
        shape is (number of types, number of sites)
    :param types: array of types, shape is (number of types,)
    :param concentrations: array of concentrations in at.%, shape is (number of types,)
    :param enthalpy_per_atom: enthalpy per atom of equilibrated configuration
    :param beta_vals: array of beta values to evaluate characteristics at
        shape is (number of temperatures,)
    :return: 3-tuple of vacancy characteristics
        each member of tuple has shape (number of temperatures,)
    """

    assert occupying_energies.shape == occupying_volumes.shape
    assert (
        vacant_energies.shape == vacant_volumes.shape == (occupying_energies.shape[1],)
    )
    assert types.shape == concentrations.shape == (occupying_energies.shape[0],)

    # calculate chemical potentials
    chemical_potentials = get_chemical_potentials(
        types, occupying_energies, concentrations, enthalpy_per_atom
    )

    # calculate local formation energies and volumes as arrays of shape (num_types, num_sites)
    energetics_data = get_formation_array(
        occupying_energies, vacant_energies, chemical_potentials
    )
    volumetrics_data = get_formation_array(occupying_volumes, vacant_volumes)

    # initialize formation calculator object, calculate concentration,
    # formation energies, and formation volumes as a function of temperature
    formation_calculator = FormationCalculator(energetics_data, volumetrics_data)
    vacancy_concentration = formation_calculator.concentration_vectorized(beta_vals)
    formation_energies = formation_calculator.get_formation_energy(beta_vals)
    formation_volumes = formation_calculator.get_formation_volume(beta_vals)

    # return tuple of formation characteristics
    return vacancy_concentration, formation_energies, formation_volumes
