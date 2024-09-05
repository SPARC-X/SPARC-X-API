"""A minimal example combining active learning library like Finetuna with SPARC

usage
First download the checkpoint from the url https://dl.fbaipublicfiles.com/opencatalystproject/models/2021_08/s2ef/gemnet_t_direct_h512_all.pt

python example_finetuna_minimal.py
"""
import argparse
import os
from pathlib import Path

import ase
import torch
import yaml
from ase.build import molecule
from ase.cluster.cubic import FaceCenteredCubic
from ase.constraints import FixAtoms
from ase.io.trajectory import Trajectory
from ase.optimize import BFGS
from finetuna.ml_potentials.finetuner_ensemble_calc import FinetunerEnsembleCalc
from finetuna.online_learner.online_learner import OnlineLearner

from sparc.calculator import SPARC

cpu = not torch.cuda.is_available()
curdir = Path(__file__).parent
config_file = curdir / "ft_config_gemnet_gpu.yml"
with open(config_file, "r") as fd:
    configs = yaml.load(fd, Loader=yaml.FullLoader)

checkpoint = os.environ.get("CHECKPOINT_PATH", None)
if checkpoint is None:
    # Use default (relative path)
    checkpoint = curdir / configs["ocp"]["checkpoint_path_list"][0]
checkpoint = Path(checkpoint)

if not checkpoint.is_file():
    raise FileNotFoundError("Cannot found the model checkpoint file!")

finetuner = configs["finetuner"]
finetuner[0].update(cpu=cpu)
learner = configs["learner"]

ml_potential = FinetunerEnsembleCalc(
    checkpoint_paths=[checkpoint],
    mlp_params=finetuner,
)


# init_molecule = molecule("H2O", pbc=False, cell=[8, 8, 8])
# init_molecule.center()
# init_molecule.rattle()


surfaces = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]
layers = [1, 2, 1]
lc = 3.61000

# init_atoms = molecule("CH4", pbc=True, cell=[8, 8, 8])
# init_atoms.constraints = [FixAtoms([0])]
# init_atoms.center()
# init_atoms.rattle(0.05)

init_atoms = FaceCenteredCubic("Cu", surfaces, layers, latticeconstant=lc)
init_atoms.cell = [12, 12, 12]
init_atoms.center()
init_atoms.pbc = False
init_atoms.rattle(0.05)


sparc_params = {"xc": "pbe", "h": 0.13}
# with SPARC(directory="pure_BFGS", **sparc_params) as calc:
#     atoms = init_atoms.copy()
#     atoms.calc = calc
#     dyn = BFGS(atoms, maxstep=0.2, trajectory="pure_bfgs.traj")
#     dyn.run(fmax=0.03)


with SPARC(directory="online_coldstart", **sparc_params) as parent_calc:
    atoms = init_atoms.copy()
    onlinecalc = OnlineLearner(learner, [], ml_potential, parent_calc)
    atoms.calc = onlinecalc
    dyn = BFGS(atoms, maxstep=0.2)
    dyn.run(fmax=0.03)
