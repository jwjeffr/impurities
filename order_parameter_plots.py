#!/usr/bin/env python

from itertools import combinations_with_replacement
import json
from logging.config import valid_ident

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import ovito
from cowley_sro_parameters import sro_modifier
from matplotlib.gridspec import GridSpec
from scoreBasedDenoising import ScoreBasedDenoising

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
    with open("config.json", "r") as file:
        config = json.load(file)

    # need to convert type map keys to integers for SRO modifier to work correctly
    type_maps = config["Type Maps"]
    for system in config["Systems"]:
        type_maps[system] = {int(key): config["Atom Abbreviations"][val] for key, val in type_maps[system].items()}

    num_systems = len(config["Systems"])

    max_num_types = max([len(map_) for map_ in config["Type Maps"].values()])

    frames = np.arange(config["Number of Frames"])
    timestep = np.zeros(config["Number of Frames"], dtype=int)
    sro_params = np.zeros(
        (num_systems, config["Number of Frames"], max_num_types, max_num_types)
    )
    frobenius_norms = np.zeros((config["Number of Frames"], 2))

    for system_index, system in enumerate(config["Systems"]):
        pipeline = ovito.io.import_file(f"mc_data/{system}/mc.dump")
        structure = config["Structure"][system]
        pipeline.modifiers.append(ScoreBasedDenoising(structure=structure))
        bonds_modifier = nearest_neighbor_topology_modifier(
            config["Number of Nearest Neighbors"][structure]
        )
        pipeline.modifiers.append(bonds_modifier)

        type_map = type_maps[system]
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
            t = type_maps["cantor"][i + 1], type_maps["cantor"][j + 1]
            ax = plt.subplot(gs[i, j])
            plot_window(
                ax,
                i,
                j,
                timestep,
                max_num_types - 1,
                t,
                sro_params[0, :, j, i],
                color=config["System Colors"]["cantor"],
                linestyle=config["System Line Styles"]["cantor"],
            )

    for i in range(2):
        for j in range(i + 1):
            t = type_maps["FeAl"][i + 1], type_maps["FeAl"][j + 1]
            ax = plt.subplot(gs[i + 6, j])
            plot_window(
                ax,
                i,
                j,
                timestep,
                1,
                t,
                sro_params[1, :, j, i],
                color=config["System Colors"]["FeAl"],
                linestyle=config["System Line Styles"]["FeAl"],
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


if __name__ == "__main__":
    main()
