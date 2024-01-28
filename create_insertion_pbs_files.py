import numpy as np

from constants import SYSTEMS


def write_file(t: int, tag: str) -> list:
    lines = [
        "#!/bin/bash",
        "#PBS -S /bin/bash",
        f"#PBS -N ins{t}{tag}",
        "#PBS -l select=1:ncpus=32:mpiprocs=32:mem=64gb:interconnect=any,walltime=72:00:00",
        "#PBS -j oe",
        "cd $PBS_O_WORKDIR",
        f"mpirun -np 32 ~/software/lammps/lammps-2Aug2023/build/lmp -in insertions.in -var tag {tag} -var t {t} -log logs/{tag}/insertions{t}.log",
    ]

    file_name = f"insertions_{tag}_{t}.pbs"

    with open(file_name, "w") as file:
        file.writelines([f"{line}\n" for line in lines])


def main():
    t_values = np.loadtxt("time.txt", dtype=int)

    for system in SYSTEMS:
        for t in t_values:
            write_file(t, system)


if __name__ == "__main__":
    main()
