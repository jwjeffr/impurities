import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from constants import SYSTEMS, SystemColors, SystemLineStyles


def temp_beta_conversion(x: float) -> float:
    boltzmann_constant = 8.615e-5

    return 1 / (boltzmann_constant * x)


def main():
    mpl.use("Agg")

    beta_vals = np.linspace(10.0, 30.0, 10_000)

    for system in SYSTEMS:
        chemical_potentials = np.loadtxt(f"chemical_potentials_{system}.txt")

        fugacities = np.exp(np.outer(beta_vals, chemical_potentials))
        squared_sum = np.sum(fugacities, axis=1) ** 2
        sum_squared = np.sum(fugacities**2, axis=1)

        plot_kwargs = {
            "color": getattr(SystemColors(), system.upper()),
            "linestyle": getattr(SystemLineStyles(), system.upper()),
            "label": system,
        }

        fluctuations = np.sqrt(1 - sum_squared / squared_sum)
        plt.plot(beta_vals, fluctuations, **plot_kwargs)

    plt.grid()
    plt.legend()
    plt.yscale("log")

    ax = plt.gca()

    secx = ax.secondary_xaxis(
        "top", functions=(temp_beta_conversion, temp_beta_conversion)
    )
    secx.set_xticks([temp_beta_conversion(x) for x in ax.get_xticks()])
    new_labels = [f"{x / 100:.2f}" for x in secx.get_xticks()]
    secx.set_xticklabels(new_labels)
    secx.set_xlabel("temperature ($10^2$ K)")

    plt.xlabel(r"inverse temperature ($\beta$) (eV$^{-1}$)")
    plt.ylabel(r"occupation number fluctuation ($\Delta n$)")

    plt.tight_layout()
    plt.savefig("plots/fluctuation.svg", bbox_inches="tight")


if __name__ == "__main__":
    main()
