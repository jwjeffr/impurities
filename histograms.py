#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import site_statistics
from constants import AtomColors, SYSTEMS, TYPE_MAPS, FINAL_STEP


def main():
    mpl.use("Agg")
    height_ratios = [1, 1, 1, 1, 1, 0.4, 1, 1]
    fig, axs = plt.subplots(
        sharex="col",
        ncols=2,
        nrows=8,
        sharey=True,
        height_ratios=height_ratios,
        figsize=(5, 5),
    )

    for system in SYSTEMS:
        type_dict = TYPE_MAPS[system]
        num_types = len(type_dict)

        types = np.arange(num_types, dtype=int)
        concentrations = np.ones(types.shape) / len(types)

        occupying_energies = np.vstack(
            [
                np.loadtxt(
                    f"energetics_data/{system}/occupying{t + 1:.0f}_{FINAL_STEP:.0f}.txt"
                )
                for t in types
            ]
        )
        vacant_energies = np.loadtxt(
            f"energetics_data/{system}/vacant_{FINAL_STEP:.0f}.txt"
        )
        enthalpy_per_atom = np.loadtxt(
            f"energetics_data/{system}/enthalpy_{FINAL_STEP:.0f}.txt"
        )

        occupying_volumes = np.vstack(
            [
                np.loadtxt(
                    f"volumetrics_data/{system}/occupying{t + 1:.0f}_{FINAL_STEP:.0f}.txt"
                )
                for t in types
            ]
        )
        vacant_volumes = np.loadtxt(
            f"volumetrics_data/{system}/vacant_{FINAL_STEP:.0f}.txt"
        )

        chemical_potentials = site_statistics.get_chemical_potentials(
            types, occupying_energies, concentrations, enthalpy_per_atom
        )
        for i, chemical_potential in enumerate(chemical_potentials):
            print(
                f"chemical potential of {type_dict[i + 1]} in {system} = {chemical_potential:.2f}"
            )
        np.savetxt(f"chemical_potentials_{system}.txt", chemical_potentials)
        formation_enthalpies = site_statistics.get_formation_array(
            occupying_energies, vacant_energies, chemical_potentials
        )
        formation_volumes = site_statistics.get_formation_array(
            occupying_volumes, vacant_volumes
        )

        plot_keyword_args = dict(zorder=6, linewidth=1, edgecolor="black")

        if system == "cantor":
            colors = [
                AtomColors.COBALT,
                AtomColors.NICKEL,
                AtomColors.CHROMIUM,
                AtomColors.IRON,
                AtomColors.MANGANESE,
            ]
        elif system == "FeAl":
            colors = [AtomColors.IRON, AtomColors.ALUMINUM]

        for i in np.arange(num_types):
            axis_index = i
            if system == "FeAl":
                axis_index += 6

            for ax in axs[axis_index]:
                ax.grid()
                ax.set_ylim([None, 550])

            enthalpy_ax, volume_ax = axs[axis_index]
            enthalpy_ax.hist(
                formation_enthalpies[i, :], **plot_keyword_args, color=(colors[i], 1.0)
            )
            volume_ax.hist(
                formation_volumes[i, :], **plot_keyword_args, color=(colors[i], 1.0)
            )
            volume_ax.text(
                0.0, 250, r"$\alpha = $" + type_dict[i + 1], va="center", ha="left"
            )

    axs[5, 0].axis("off")
    axs[5, 1].axis("off")
    axs[-1, 0].set_xlabel(
        "local formation enthalpy\n" + r"$\mathcal{H}_\sigma^{(\alpha)}$ (eV)"
    )
    axs[-1, 1].set_xlabel(
        "local formation volume\n" + r"$v_\sigma - V_\sigma^{(\alpha)}$ ($\AA^3$)"
    )

    fig.text(0.0, 0.5, "counts", ha="left", va="center", rotation="vertical")
    fig.tight_layout(w_pad=2)
    fig.savefig("plots/distribution.svg")


if __name__ == "__main__":
    main()
