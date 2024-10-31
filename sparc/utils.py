"""Utilities that are loosely related to core sparc functionalities
"""
import _thread
import io
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
from contextlib import contextmanager
from pathlib import Path
from typing import List, Optional, Union
from warnings import warn

import numpy as np
import psutil

from .api import SparcAPI
from .docparser import SparcDocParser


def deprecated(message):
    def decorator(func):
        def new_func(*args, **kwargs):
            warn(
                "Function {} is deprecated! {}".format(func.__name__, message),
                category=DeprecationWarning,
            )
            return func(*args, **kwargs)

        return new_func

    return decorator


def compare_dict(d1, d2):
    """Helper function to compare dictionaries"""
    # Use symmetric difference to find keys which aren't shared
    # for python 2.7 compatibility
    if set(d1.keys()) ^ set(d2.keys()):
        return False

    # Check for differences in values
    for key, value in d1.items():
        if np.any(value != d2[key]):
            return False
    return True


def string2index(string: str) -> Union[int, slice, str]:
    """Convert index string to either int or slice
    This method is a copy of ase.io.formats.string2index
    """
    # A quick fix for slice
    if isinstance(string, (list, slice)):
        return string
    if ":" not in string:
        # may contain database accessor
        try:
            return int(string)
        except ValueError:
            return string
    i: List[Optional[int]] = []
    for s in string.split(":"):
        if s == "":
            i.append(None)
        else:
            i.append(int(s))
    i += (3 - len(i)) * [None]
    return slice(*i)


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


def locate_api(json_file=None, doc_path=None):
    """Find the default api in the following order
    1) User-provided json file path
    2) User-provided path to the doc
    3) If none of the above is provided, try to use SPARC_DOC_PATH
    4) Fallback to the as-shipped json api
    """
    if json_file is not None:
        api = SparcAPI(json_file)
        return api

    if doc_path is None:
        doc_path = os.environ.get("SPARC_DOC_PATH", None)

    if (doc_path is not None) and Path(doc_path).is_dir():
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmpdir = Path(tmpdir)
                tmpfile = tmpdir / "parameters.json"
                with open(tmpfile, "w") as fd:
                    fd.write(
                        SparcDocParser.json_from_directory(
                            Path(doc_path), include_subdirs=True
                        )
                    )
                api = SparcAPI(tmpfile)
            api.source["path"] = Path(doc_path).resolve().as_posix()
            api.source["type"] = "latex"
            return api
        except Exception as e:
            warn(f"Cannot load JSON schema from env {doc_path}, the error is {e}.")
            pass

    api = SparcAPI()
    return api


# Utilities taken from vasp_interactive project


class TimeoutException(Exception):
    """Simple class for timeout"""

    pass


@contextmanager
def time_limit(seconds):
    """Usage:
    try:
        with time_limit(60):
            do_something()
    except TimeoutException:
        raise
    """

    def signal_handler(signum, frame):
        raise TimeoutException("Timed out closing sparc process.")

    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)


class ProcessReturned(Exception):
    """Simple class for process that has returned"""

    pass


@contextmanager
def monitor_process(self, interval=1.0):
    """Usage:
    try:
        with monitor_process(process):
            do_something()
    except TimeoutException:
        raise
    """

    def signal_handler(signum, frame):
        raise ProcessReturned(
            f"Process {self.process.pid} has returned with exit code {self.process.poll()}!"
        )

    def check_process():
        while True:
            if self.process.poll() is not None:
                # signal.alarm(0)
                print("The process has exited")
                self.in_socket.close()
                print(self.in_socket)
                signal(signal.SIGALRM)
                raise ProcessReturned(
                    f"Process {self.process.pid} has returned with exit code {self.process.poll()}!"
                )
            time.sleep(interval)

    if self.process is None:
        raise RuntimeError("No process selected!")

    signal.signal(signal.SIGALRM, signal_handler)
    monitor = threading.Thread(target=check_process)
    monitor.start()
    try:
        yield
    finally:
        monitor.join()


def _find_mpi_process(pid, mpi_program="mpirun", sparc_program="sparc"):
    """Recursively search children processes with PID=pid and return the one
    that mpirun (or synonyms) are the main command.

    If srun is found as the process, need to use `scancel` to pause / resume the job step
    """
    allowed_names = set(["mpirun", "mpiexec", "orterun", "oshrun", "shmemrun"])
    allowed_sparc_names = set(["sparc"])
    if mpi_program:
        allowed_names.add(mpi_program)
    if sparc_program:
        allowed_sparc_names.add(sparc_program)
    try:
        process_list = [psutil.Process(pid)]
    except psutil.NoSuchProcess:
        warn(
            "Psutil cannot locate the pid. Your sparc program may have already exited."
        )
        match = {"type": None, "process": None}
        return match

    process_list.extend(process_list[0].children(recursive=True))
    mpi_candidates = []
    match = {"type": None, "process": None}
    for proc in process_list:
        name = proc.name()
        if name in ["srun"]:
            match["type"] = "slurm"
            match["process"] = _locate_slurm_step(program=sparc_program)
            break
        elif proc.name() in allowed_names:
            # are the mpi process's direct children sparc binaries?
            children = proc.children()
            if len(children) > 0:
                if children[0].name() in allowed_sparc_names:
                    mpi_candidates.append(proc)
    if len(mpi_candidates) > 1:
        warn(
            "More than 1 mpi processes are created. This may be a bug. I'll use the last one"
        )
    if len(mpi_candidates) > 0:
        match["type"] = "mpi"
        match["process"] = mpi_candidates[-1]

    return match


def _get_slurm_jobid():
    jobid = os.environ.get("SLURM_JOB_ID", None)
    if jobid is None:
        jobid = os.environ.get("SLURM_JOBID", None)
    return jobid


def _locate_slurm_step(program="sparc"):
    """If slurm job system is involved, search for the slurm step id
    that matches vasp_std (or other vasp commands)

    Steps:
    1. Look for SLURM_JOB_ID in current env
    2. Use `squeue` to locate the sparc step (latest)

    squeue
    """
    allowed_names = set(["sparc"])
    if program:
        allowed_names.add(program)
    jobid = _get_slurm_jobid()
    if jobid is None:
        # TODO: convert warn to logger
        warn(("Cannot locate the SLURM job id."))
        return None
    # Only 2 column output (jobid and jobname)
    cmds = ["squeue", "-s", "--job", str(jobid), "-o", "%.30i %.30j"]
    proc = _run_process(cmds, capture_output=True)
    output = proc.stdout.decode("utf8").split("\n")
    # print(output)
    candidates = []
    # breakpoint()
    for line in output[1:]:
        try:
            stepid, name = line.strip().split()
        except Exception:
            continue
        if any([v in name for v in allowed_names]):
            candidates.append(stepid)

    if len(candidates) > 1:
        warn("More than 1 slurm steps are found. I'll use the most recent one")
    if len(candidates) > 0:
        proc = candidates[0]
    else:
        proc = None
    return proc


def _slurm_signal(stepid, sig=signal.SIGTSTP):
    if isinstance(sig, (str,)):
        sig = str(sig)
    elif isinstance(sig, (int,)):
        sig = signal.Signals(sig).name
    else:
        sig = sig.name
    cmds = ["scancel", "-s", sig, str(stepid)]
    proc = _run_process(cmds, capture_output=True)
    output = proc.stdout.decode("utf8").split("\n")
    return


def _run_process(commands, shell=False, print_cmd=True, cwd=".", capture_output=False):
    """Wrap around subprocess.run
    Returns the process object
    """
    full_cmd = " ".join(commands)
    if print_cmd:
        print(" ".join(commands))
    if shell is False:
        proc = subprocess.run(
            commands, shell=shell, cwd=cwd, capture_output=capture_output
        )
    else:
        proc = subprocess.run(
            full_cmd, shell=shell, cwd=cwd, capture_output=capture_output
        )
    if proc.returncode == 0:
        return proc
    else:
        raise RuntimeError(f"Running {full_cmd} returned error code {proc.returncode}")
