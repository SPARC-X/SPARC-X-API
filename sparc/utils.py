"""Utilities that are loosely related to core sparc functionalities
"""
import os
import shutil
import numpy as np


def _find_default_sparc():
    """Find the default sparc by $PATH and mpi location"""
    sparc_exe = shutil.which("sparc")

    mpi_exe = shutil.which("mpirun")
    # TODO: more examples on pbs / lsf
    if mpi_exe is not None:
        try:
            num_cores = int(
                os.environ.get(
                    "OMPI_COMM_WORLD_SIZE",
                    os.environ.get(
                        "OMPI_UNIVERSE_SIZE",
                        os.environ.get("MPICH_RANK_REORDER_METHOD", ""),
                    ).split(":")[-1],
                )
            )
        except Exception:
            num_cores = 1
        return sparc_exe, mpi_exe, num_cores

    mpi_exe = shutil.which("srun")
    if mpi_exe is not None:
        # If srun is available, get the number of cores from the environment
        num_cores = int(os.environ.get("SLURM_JOB_CPUS_PER_NODE", 1))
        return sparc_exe, mpi_exe, num_cores

    return sparc_exe, None, 1


def h2gpts(h, cell_cv, idiv=4):
    """Convert a h-parameter (Angstrom) to gpts"""
    cell_cv = np.array(cell_cv)
    cell_lengths = np.linalg.norm(cell_cv, axis=1)
    grid = np.ceil(cell_lengths / h)
    grid = np.maximum(idiv, grid)
    return [int(a) for a in grid]


def cprint(content, color=None, bold=False, underline=False, **kwargs):
    """Color print wrapper for ansi terminal.
    Only a few color names are provided
    """
    ansi_color = dict(
        HEADER="\033[95m",
        COMMENT="\033[90m",
        OKBLUE="\033[94m",
        OKGREEN="\033[92m",
        OKCYAN="\033[96m",
        WARNING="\033[93m",
        FAIL="\033[91m",
        ENDC="\033[0m",
    )

    style_codes = {"BOLD": "\033[1m", "UNDERLINE": "\033[4m"}

    if color is None:
        output = content
    elif color.upper() in ansi_color.keys() and color.upper() != "ENDC":
        output = ansi_color[color.upper()] + content + ansi_color["ENDC"]
    else:
        raise ValueError(
            f"Unknown ANSI color name. Allowed values are {list(ansi_color.keys())}"
        )

    if bold:
        output = style_codes["BOLD"] + output + ansi_color["ENDC"]

    if underline:
        output = style_codes["UNDERLINE"] + output + ansi_color["ENDC"]

    print(output, **kwargs)
    return
