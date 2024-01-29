#!/usr/bin/python

"""
script for creating formation plots at the final MC-MD step
"""

import json

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
            "label": system,
        }

        axs[0].plot(beta_vals, vacancy_concentration, **plot_kwargs)
        axs[1].plot(beta_vals, formation_energies, **plot_kwargs)
        axs[2].plot(beta_vals, formation_volumes, **plot_kwargs)

    axs[0].set_yscale("log")
    axs[0].grid()
    axs[0].set_ylabel(r"concentration $x_V$ (at. %)")
    axs[0].legend()

    axs[1].grid()
    axs[1].set_ylabel(r"formation energy $E_{form}$ (eV)")

    axs[2].grid()
    axs[2].set_ylabel(r"formation volume $\Omega_{form}$ ($\AA^3$)")

    axs[-1].set_xlabel(r"inverse temperature ($\beta$) (eV$^{-1}$)")

    secx = axs[0].secondary_xaxis(
        "top", functions=(temp_beta_conversion, temp_beta_conversion)
    )
    secx.set_xticks([temp_beta_conversion(x) for x in axs[0].get_xticks()])
    new_labels = [f"{x / 100:.2f}" for x in secx.get_xticks()]
    secx.set_xticklabels(new_labels)
    secx.set_xlabel("temperature ($10^2$ K)")

    fig.tight_layout()
    fig.savefig("plots/formations.svg", bbox_inches="tight")


if __name__ == "__main__":
    main()
