#!/usr/bin/python

"""
script for creating formation plots at the final MC-MD step
"""

import json
import re

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from site_statistics import get_vacancy_characteristics


def temp_beta_conversion(x: float) -> float:
    """
    temperature <-> beta conversion
    """

    boltzmann_constant = 8.615e-5

    return 1 / (boltzmann_constant * x)


def main():
    """
    create plot
    """

    mpl.use("Agg")

    with open("config.json", "r", encoding="utf8") as file:
        config = json.load(file)

    fig, axs = plt.subplots(sharex=True, nrows=3, figsize=(6, 8))

    for system in config["Systems"]:
        type_dict = config["Type Maps"][system]
        num_types = len(type_dict)

        types = np.arange(num_types, dtype=int)
        concentrations = np.ones(types.shape) / len(types)

        occupying_energies = np.vstack(
            [
                np.loadtxt(
                    f"energetics_data/{system}/occupying{t + 1:.0f}_{config['Final Step']:.0f}.txt"
                )
                for t in types
            ]
        )
        vacant_energies = np.loadtxt(
            f"energetics_data/{system}/vacant_{config['Final Step']:.0f}.txt"
        )
        enthalpy_per_atom = np.loadtxt(
            f"energetics_data/{system}/enthalpy_{config['Final Step']:.0f}.txt"
        )

        occupying_volumes = np.vstack(
            [
                np.loadtxt(
                    f"volumetrics_data/{system}/occupying{t + 1:.0f}_{config['Final Step']:.0f}.txt"
                )
                for t in types
            ]
        )
        vacant_volumes = np.loadtxt(
            f"volumetrics_data/{system}/vacant_{config['Final Step']:.0f}.txt"
        )

        beta_vals = np.linspace(10.0, 30.0, 10_000)
        (
            vacancy_concentration,
            formation_energies,
            formation_volumes,
        ) = get_vacancy_characteristics(
            vacant_energies,
            occupying_energies,
            vacant_volumes,
            occupying_volumes,
            types,
            concentrations,
            enthalpy_per_atom,
            beta_vals,
        )

        plot_kwargs = {
            "color": config["System Colors"][system],
            "linestyle": config["System Line Styles"][system],
            "label": config["System Labels"][system],
        }

        axs[0].plot(beta_vals, vacancy_concentration, **plot_kwargs)
        axs[1].plot(beta_vals, formation_energies, **plot_kwargs)
        axs[2].plot(beta_vals, formation_volumes, **plot_kwargs)

        # compute two-state system parameters
        # kind of hacky but two-state isn't purpose of work, so not ingrained in any libraries

        num_atoms, dx, dy, dz = None, None, None, None

        with open(
            f"mc_data/{system}/mc_relaxed_{config['Final Step']:.0f}.dat", "r"
        ) as file:
            for line in file:
                if num_atoms and dx and dy and dz:
                    continue
                if (match := re.match(r"(\d+) atoms", line)) and not num_atoms:
                    num_atoms = int(match.group(1))
                elif (
                    match := re.match(r"(-?\d+\.\d+?) (-?\d+\.\d+?) xlo xhi", line)
                ) and not dx:
                    xlo, xhi = match.groups()
                    dx = float(xhi) - float(xlo)
                elif (
                    match := re.match(r"(-?\d+\.\d+?) (-?\d+\.\d+?) ylo yhi", line)
                ) and not dy:
                    ylo, yhi = match.groups()
                    dy = float(yhi) - float(ylo)
                elif (
                    match := re.match(r"(-?\d+\.\d+?) (-?\d+\.\d+?) zlo zhi", line)
                ) and not dz:
                    zlo, zhi = match.groups()
                    dz = float(zhi) - float(zlo)

        if not (num_atoms and dx and dy and dz):
            raise ValueError

        vacant_energies = np.loadtxt(
            f"energetics_data/{system}/vacant_{config['Final Step']:.0f}.txt"
        )
        vacant_volumes = np.loadtxt(
            f"volumetrics_data/{system}/vacant_{config['Final Step']:.0f}.txt"
        )
        enthalpy_per_atom = np.loadtxt(
            f"energetics_data/{system}/enthalpy_{config['Final Step']:.0f}.txt"
        )

        # should be skipping a different number of rows for numpy <= 1.22
        # https://numpy.org/devdocs/release/1.23.0-notes.html
        labels, occupying_types = np.loadtxt(
            f"mc_data/{system}/mc_relaxed_{config['Final Step']:.0f}.dat",
            skiprows=13 + num_types,
            max_rows=num_atoms,
            usecols=[0, 1],
            dtype=int,
        ).T
        chemical_potentials = np.loadtxt(f"chemical_potentials_{system}.txt")
        chemical_potentials = chemical_potentials[
            occupying_types[np.argsort(labels)] - 1
        ]

        reference_energy = num_atoms * enthalpy_per_atom
        reference_volume = dx * dy * dz

        formation_energies = vacant_energies - reference_energy + chemical_potentials
        formation_volumes = vacant_volumes - reference_volume

        vacant_probabilities = 1.0 / (
            1.0 + np.exp(np.outer(beta_vals, formation_energies))
        )
        concentration = np.mean(vacant_probabilities, axis=1)

        numerator = formation_energies * np.exp(np.outer(beta_vals, formation_energies))
        denominator = (1.0 + np.exp(np.outer(beta_vals, formation_energies))) ** 2
        formation_energy = (
            1.0 / concentration * np.mean(numerator / denominator, axis=1)
        )

        numerator = formation_volumes * np.exp(np.outer(beta_vals, formation_energies))
        formation_volume = (
            1.0 / concentration * np.mean(numerator / denominator, axis=1)
        )

        plot_kwargs["label"] = f"{plot_kwargs['label']} (two-state)"
        plot_kwargs["linestyle"] = ":"
        axs[0].plot(beta_vals, concentration, **plot_kwargs)
        axs[1].plot(beta_vals, formation_energy, **plot_kwargs)
        axs[2].plot(beta_vals, formation_volume, **plot_kwargs)

    axs[0].set_yscale("log")
    axs[0].grid()
    axs[0].set_ylabel(r"concentration $x_V$ (at. fraction)")
    axs[0].legend()

    axs[1].grid()
    axs[1].set_ylabel(r"formation energy $E_{form}$ (eV)")

    axs[2].grid()
    axs[2].set_ylabel(r"formation volume $\Omega_{form}$ ($\AA^3$)")

    axs[-1].set_xlabel(r"inverse temperature ($\beta$) (eV$^{-1}$)")

    secx = axs[0].secondary_xaxis(
        "top", functions=(temp_beta_conversion, temp_beta_conversion)
    )
    temperature_spacing = config["Temperature Spacing"]
    min_temperature = temperature_spacing * round(
        temp_beta_conversion(max(axs[0].get_xticks())) / temperature_spacing
    )
    max_temperature = temperature_spacing * round(
        temp_beta_conversion(min(axs[0].get_xticks())) / temperature_spacing
    )
    secx.set_xticks(
        np.arange(
            min_temperature,
            max_temperature + temperature_spacing,
            step=temperature_spacing,
        )
    )
    new_labels = [f"{x / 100:.0f}" for x in secx.get_xticks()]
    secx.set_xticklabels(new_labels)
    secx.set_xlabel("temperature ($10^2$ K)")

    fig.tight_layout()
    fig.savefig("plots/formations.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
