import os
import shutil
import numpy as np


def _find_default_sparc():
    """Find the default sparc by $PATH and mpi location"""
    sparc_exe = shutil.which("sparc")

    mpi_exe = shutil.which("mpirun")
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
