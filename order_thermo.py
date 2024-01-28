#!/usr/bin/env python

import json

import ovito
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from cowley_sro_parameters import sro_modifier
from scoreBasedDenoising import ScoreBasedDenoising

from modifiers import nearest_neighbor_topology_modifier
from site_statistics import get_vacancy_characteristics


def temp_beta_conversion(x: float) -> float:
    boltzmann_constant = 8.615e-5

    return 1 / (boltzmann_constant * x)


def main():
    mpl.use("Agg")
    with open("config.json", "r") as file:
        config = json.load(file)

    temperatures = np.arange(200, 1000 + 200, step=200)

    multipanel_fig, axs = plt.subplots(
        nrows=3, ncols=len(config["Systems"]), sharex="col"
    )

    color_map = mpl.colormaps["inferno"]
    normalized_temperatures = (temperatures - min(temperatures)) / (
        max(temperatures) - min(temperatures)
    )
    color_map = color_map(normalized_temperatures)

    for system_index, system in enumerate(config["Systems"]):
        type_map = config["Type Maps"][system]

        type_map = {
            int(key): config["Atom Abbreviations"][val] for key, val in type_map.items()
        }

        num_types = len(type_map)
        types = np.arange(num_types)
        concentrations = np.ones_like(types) / num_types

        order_parameters = np.zeros(config["Number of Frames"])
        timesteps = np.zeros(config["Number of Frames"], dtype=int)

        mc_file_name = f"mc_data/{system}/mc.dump"
        pipeline = ovito.io.import_file(mc_file_name)

        structure = config["Structure"][system]
        attribute = f"sro_{config['Dominant Order Parameter'][system]}"
        pipeline.modifiers.append(ScoreBasedDenoising(structure=structure))

        bonds_modifier = nearest_neighbor_topology_modifier(
            config["Number of Nearest Neighbors"][structure]
        )
        pipeline.modifiers.append(bonds_modifier)

        modifier = sro_modifier(type_map=type_map)
        pipeline.modifiers.append(modifier)

        for frame in np.arange(config["Number of Frames"]):
            data = pipeline.compute(frame)
            timesteps[frame] = data.attributes["Timestep"]
            order_parameters[frame] = data.attributes[attribute]

        for temperature, color in zip(temperatures, color_map):
            vacancy_concentrations = np.zeros(config["Number of Frames"])
            formation_energies = np.zeros(config["Number of Frames"])
            formation_volumes = np.zeros(config["Number of Frames"])

            beta = temp_beta_conversion(temperature)
            plot_kwargs = {
                "edgecolor": "black",
                "facecolor": color,
                "alpha": 0.7,
                "zorder": 6,
            }

            for frame, step in enumerate(timesteps):
                occupying_energies = np.vstack(
                    [
                        np.loadtxt(
                            f"energetics_data/{system}/occupying{t + 1:.0f}_{step:.0f}.txt"
                        )
                        for t in types
                    ]
                )
                vacant_energies = np.loadtxt(
                    f"energetics_data/{system}/vacant_{step:.0f}.txt"
                )
                enthalpy_per_atom = np.loadtxt(
                    f"energetics_data/{system}/enthalpy_{step:.0f}.txt"
                )

                occupying_volumes = np.vstack(
                    [
                        np.loadtxt(
                            f"volumetrics_data/{system}/occupying{t + 1:.0f}_{step:.0f}.txt"
                        )
                        for t in types
                    ]
                )
                vacant_volumes = np.loadtxt(
                    f"volumetrics_data/{system}/vacant_{step:.0f}.txt"
                )

                (
                    vacancy_concentration,
                    formation_energy,
                    formation_volume,
                ) = get_vacancy_characteristics(
                    vacant_energies,
                    occupying_energies,
                    vacant_volumes,
                    occupying_volumes,
                    types,
                    concentrations,
                    enthalpy_per_atom,
                    beta,
                )

                vacancy_concentrations[frame] = vacancy_concentration
                formation_energies[frame] = formation_energy
                formation_volumes[frame] = formation_volume

            axs[0, system_index].scatter(
                order_parameters,
                vacancy_concentrations,
                label=f"{temperature:.0f}",
                **plot_kwargs,
            )
            axs[1, system_index].scatter(
                order_parameters, formation_energies, **plot_kwargs
            )
            axs[2, system_index].scatter(
                order_parameters, formation_volumes, **plot_kwargs
            )

        axs[-1, system_index].set_xlabel(r"$\chi_{" + attribute + r"}$")
        axs[0, system_index].set_title(system)

    for i, ax_list in enumerate(axs):
        for ax in ax_list:
            if i == 0:
                ax.set_yscale("log")
            ax.grid()

    axs[0, -1].legend()
    handles, labels = axs[0, -1].get_legend_handles_labels()
    axs[0, -1].legend(
        handles[::-1], labels[::-1], title="temperature (K)", bbox_to_anchor=(1.1, 1.1)
    )

    axs[0, 0].set_ylabel(r"$x_V$ (at. %)")
    axs[1, 0].set_ylabel(r"$E_{form}$ (eV)")
    axs[2, 0].set_ylabel(r"$\Omega_{form}$ ($\AA^3$)")

    multipanel_fig.tight_layout()
    multipanel_fig.savefig("plots/order_thermo.svg", bbox_inches="tight")


if __name__ == "__main__":
    main()
