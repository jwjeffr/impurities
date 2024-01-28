from itertools import combinations_with_replacement

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import ovito
from cowley_sro_parameters import sro_modifier
from matplotlib.gridspec import GridSpec
from scoreBasedDenoising import ScoreBasedDenoising

from constants import (
    SYSTEMS,
    NUM_FRAMES,
    TYPE_MAPS,
    NEAREST_NEIGHBORS,
    SystemColors,
    SystemLineStyles,
)
from modifiers import nearest_neighbor_topology_modifier


def plot_window(ax, i, j, timestep, bottom: int, t: tuple, params, color, linestyle):
    ax.set_xlim([10 ** (0.5), 10 ** (6.5)])
    ax.set_ylim([-1.1, 1.1])
    ax.set_xscale("log")
    ax.set_xticks([1e1, 1e2, 1e3, 1e4, 1e5, 1e6])
    ax.set_yticks([-1, -0.5, 0, 0.5, 1])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.grid()
    # along the left side of the grid, label y-axis ticks
    if j == 0:
        ax.set_yticklabels([-1, "", 0, "", 1])

    # along the bottom side of the grid, label x-axis ticks
    if i == bottom:
        ax.set_xticklabels(["$10^1$", "", "", "", "", "$10^6$"])

    # along the diagonal, label atom type names
    if i == j:
        ax.text(10 ** (3.5), 1.3, f"α = {t[1]}", ha="center", va="bottom")
        ax.text(10 ** (7.0), 0.0, f"α' = {t[0]}", ha="left", va="center")
    ax.plot(timestep, params, color=color, linestyle=linestyle)


def main():
    mpl.use("Agg")
    num_systems = len(SYSTEMS)

    max_num_types = max([len(map_) for map_ in TYPE_MAPS.values()])

    frames = np.arange(NUM_FRAMES)
    timestep = np.zeros(NUM_FRAMES, dtype=int)
    sro_params = np.zeros((num_systems, NUM_FRAMES, max_num_types, max_num_types))
    frobenius_norms = np.zeros((NUM_FRAMES, 2))

    for system_index, system in enumerate(SYSTEMS):
        pipeline = ovito.io.import_file(f"mc_data/{system}/mc.dump")
        if system == "cantor":
            structure = "FCC"
        if system == "FeAl":
            structure = "BCC"
        pipeline.modifiers.append(ScoreBasedDenoising(structure=structure))
        bonds_modifier = nearest_neighbor_topology_modifier(
            NEAREST_NEIGHBORS[structure]
        )
        pipeline.modifiers.append(bonds_modifier)

        type_map = TYPE_MAPS[system]
        inverse_type_map = {val: key for key, val in type_map.items()}
        pairs = list(combinations_with_replacement(type_map.values(), 2))

        # create SRO modifier which calculates all SRO's and the Frobenius-normed SRO matrix at each timestep
        modifier = sro_modifier(type_map=type_map)
        pipeline.modifiers.append(modifier)

        for frame in frames:
            data = pipeline.compute(frame)
            timestep[frame] = data.attributes["Timestep"]
            for e1, e2 in pairs:
                i, j = inverse_type_map[e1], inverse_type_map[e2]
                sro_params[system_index, frame, i - 1, j - 1] = data.attributes[
                    f"sro_{e1}{e2}"
                ]
            frobenius_norms[frame, system_index] = data.attributes[
                "frobenius_norm_sro"
            ] / len(type_map)

    np.savetxt("time.txt", timestep, fmt="%d")

    height_ratios = [1, 1, 1, 1, 1, 0.4, 1, 1]

    fig = plt.figure(figsize=(5, 8))
    gs = GridSpec(
        max_num_types + 3,
        max_num_types,
        figure=fig,
        wspace=0.3,
        hspace=0.3,
        height_ratios=height_ratios,
    )

    for i in range(max_num_types):
        for j in range(i + 1):
            t = TYPE_MAPS["cantor"][i + 1], TYPE_MAPS["cantor"][j + 1]
            ax = plt.subplot(gs[i, j])
            plot_window(
                ax,
                i,
                j,
                timestep,
                max_num_types - 1,
                t,
                sro_params[0, :, j, i],
                color=SystemColors.CANTOR,
                linestyle=SystemLineStyles.CANTOR,
            )

    for i in range(2):
        for j in range(i + 1):
            t = TYPE_MAPS["FeAl"][i + 1], TYPE_MAPS["FeAl"][j + 1]
            ax = plt.subplot(gs[i + 6, j])
            plot_window(
                ax,
                i,
                j,
                timestep,
                1,
                t,
                sro_params[1, :, j, i],
                color=SystemColors.FEAL,
                linestyle=SystemLineStyles.FEAL,
            )

    fig.text(0.5, 0.05, "MC-MD time (fs)", ha="center", va="bottom")
    fig.text(
        0.025,
        0.5,
        r"α-α' Cowley SRO parameter ($\chi_{\alpha\alpha'}$)",
        ha="left",
        va="center",
        rotation="vertical",
    )
    fig.savefig("plots/sro_params.svg", bbox_inches="tight")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for system_index, system in enumerate(SYSTEMS):
        plot_kwargs = {
            "color": getattr(SystemColors(), system.upper()),
            "linestyle": getattr(SystemLineStyles(), system.upper()),
            "label": system,
        }
        ax.plot(timestep, frobenius_norms[:, system_index], **plot_kwargs)

    ax.legend()
    ax.grid()
    ax.set_xscale("log")
    ax.set_xlabel("MC-MD time (fs)")
    ax.set_ylabel("normalized Frobenius norm of SRO matrix (Q)")
    fig.savefig("plots/frobenius.svg")


if __name__ == "__main__":
    main()
