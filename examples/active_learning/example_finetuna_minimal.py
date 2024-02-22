"""A minimal example combining active learning library like Finetuna with SPARC

usage
First download the checkpoint from the url https://dl.fbaipublicfiles.com/opencatalystproject/models/2021_08/s2ef/gemnet_t_direct_h512_all.pt

python example_finetuna_minimal.py
"""
import torch
import os
import yaml
from pathlib import Path
from ase.io.trajectory import Trajectory
from ase.optimize import BFGS
from finetuna.ml_potentials.finetuner_ensemble_calc import FinetunerEnsembleCalc
from finetuna.online_learner.online_learner import OnlineLearner
import argparse
from sparc.calculator import SPARC

from ase.build import molecule

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


init_molecule = molecule("H2O", pbc=False, cell=[8, 8, 8])

sparc_params = {"xc": "pbe", "h": 0.22}
with SPARC(**sparc_params) as parent_calc:
    onlinecalc = OnlineLearner(learner, [], ml_potential, parent_calc)
    init_molecule.calc = onlinecalc
    dyn = BFGS(init_molecule,
               maxstep=0.2)
    dyn.run(fmax=0.03)
