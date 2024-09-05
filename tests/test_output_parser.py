import os
import tempfile
from pathlib import Path

import numpy as np
import pytest
from ase.units import Bohr, Hartree, eV, kB

curdir = Path(__file__).parent
test_output_dir = curdir / "outputs"


def test_output_date_parser():
    from sparc.sparc_parsers.out import _read_sparc_version

    header1 = """***************************************************************************
*                       SPARC (version Feb 03, 2023)                      *
*   Copyright (c) 2020 Material Physics & Mechanics Group, Georgia Tech   *
*           Distributed under GNU General Public License 3 (GPL)          *
*                   Start time: Sun Feb  5 13:39:04 2023                  *
***************************************************************************
                           Input parameters
***************************************************************************"""
    assert _read_sparc_version(header1) == "2023.02.03"
    header2 = """***************************************************************************
*                       SPARC (version June 24, 2023)                      *
*   Copyright (c) 2020 Material Physics & Mechanics Group, Georgia Tech   *
*           Distributed under GNU General Public License 3 (GPL)          *
*                   Start time: Sun Feb  5 13:39:04 2023                  *
***************************************************************************
                           Input parameters
***************************************************************************"""
    assert _read_sparc_version(header2) == "2023.06.24"


def test_output_parser():
    from sparc.common import repo_dir
    from sparc.sparc_parsers.out import _read_out

    data_dict = _read_out(
        test_output_dir
        / "AlSi_primitive_quick_relax.sparc"
        / "AlSi_primitive_quick_relax.out"
    )
    assert "out" in data_dict
    out_dict = data_dict["out"]
    assert out_dict["sparc_version"] == "2023.02.03"
    # Currently the fields in run_info are strings
    assert int(out_dict["run_info"]["number of processors"]) == 48
    assert out_dict["parameters"]["NP_SPIN_PARAL"] == 1
    assert out_dict["parameters"]["LATVEC"].shape == (3, 3)
    ionic_steps = out_dict["ionic_steps"]
    assert len(ionic_steps) == 5
    for step in ionic_steps:
        assert "convergence" in step
        assert "total free energy" in step


def test_incomplete_output():
    """Test if output is incomplete, e.g. incomplete runs"""
    from sparc.sparc_parsers.out import _read_out

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        out_file = tmpdir / "test.out"
        content = """***************************************************************************
*                       SPARC (version Oct 31, 2023)                      *
*   Copyright (c) 2020 Material Physics & Mechanics Group, Georgia Tech   *
*           Distributed under GNU General Public License 3 (GPL)          *
*                   Start time: Fri Jan 19 08:30:23 2024                  *
***************************************************************************
                           Input parameters
***************************************************************************
LATVEC_SCALE: 1 1 1
LATVEC:
22.676713510043140 0.000000000000000 0.000000000000000
0.000000000000000 25.561338867158440 0.000000000000000
0.000000000000000 0.000000000000000 23.803574206414829
FD_GRID: 43 49 45
FD_ORDER: 12
BC: P P P
KPOINT_GRID: 1 1 1
KPOINT_SHIFT: 0 0 0
SPIN_TYP: 0
ELEC_TEMP_TYPE: Fermi-Dirac
SMEARING: 0.000950043469
EXCHANGE_CORRELATION: GGA_PBE
NSTATES: 9
CHEB_DEGREE: 17
CHEFSI_BOUND_FLAG: 0
CALC_STRESS: 1
MAXIT_SCF: 100
MINIT_SCF: 2
MAXIT_POISSON: 3000
TOL_SCF: 1.00E-03
POISSON_SOLVER: AAR
TOL_POISSON: 1.00E-05
TOL_LANCZOS: 1.00E-02
TOL_PSEUDOCHARGE: 1.00E-06
MIXING_VARIABLE: density
MIXING_PRECOND: kerker
TOL_PRECOND: 2.77E-04
PRECOND_KERKER_KTF: 1
PRECOND_KERKER_THRESH: 0
MIXING_PARAMETER: 1
MIXING_HISTORY: 7
PULAY_FREQUENCY: 1
PULAY_RESTART: 0
REFERENCE_CUTOFF: 0.5
RHO_TRIGGER: 4
NUM_CHEFSI: 1
FIX_RAND: 0
VERBOSITY: 1
PRINT_FORCES: 1
PRINT_ATOMS: 1
PRINT_EIGEN: 0
PRINT_DENSITY: 0
PRINT_ENERGY_DENSITY: 0
OUTPUT_FILE: SPARC
***************************************************************************
                              Socket Mode
***************************************************************************
SOCKET_HOST: localhost
SOCKET_PORT: 12345
SOCKET_INET: 1
SOCKET_MAX_NITER: 10000
***************************************************************************
                                Cell
***************************************************************************
Lattice vectors (Bohr):
22.676713510043140 0.000000000000000 0.000000000000000
0.000000000000000 25.561338867158440 0.000000000000000
0.000000000000000 0.000000000000000 23.803574206414829
Volume: 1.3797674149E+04 (Bohr^3)
Density: 1.3056802042E-03 (amu/Bohr^3), 1.4631286628E-02 (g/cc)
***************************************************************************
                           Parallelization
***************************************************************************
NP_SPIN_PARAL: 1
NP_KPOINT_PARAL: 1
NP_BAND_PARAL: 1
NP_DOMAIN_PARAL: 1 1 1
NP_DOMAIN_PHI_PARAL: 1 1 1
EIG_SERIAL_MAXNS: 1500
***************************************************************************
                             Initialization
***************************************************************************
Number of processors               :  1
Mesh spacing in x-direction        :  0.527365 (Bohr)
Mesh spacing in y-direction        :  0.52166 (Bohr)
Mesh spacing in z-direction        :  0.528968 (Bohr)
Number of symmetry adapted k-points:  1
Output printed to                  :  SPARC.out
Total number of atom types         :  2
Total number of atoms              :  3
Total number of electrons          :  8
Atom type 1  (valence electrons)   :  H 1
Pseudopotential                    :  01_H_1_1.0_1.0_pbe_v1.0.psp8
Atomic mass                        :  1.007975
Pseudocharge radii of atom type 1  :  4.22 4.17 4.23 (x, y, z dir)
Number of atoms of type 1          :  2
Atom type 2  (valence electrons)   :  O 6
Pseudopotential                    :  08_O_6_1.2_1.4_pbe_n_v1.0.psp8
Atomic mass                        :  15.9994
Pseudocharge radii of atom type 2  :  7.38 7.30 7.41 (x, y, z dir)
Number of atoms of type 2          :  1
Estimated total memory usage       :  47.35 MB
Estimated memory per processor     :  47.35 MB
===================================================================
                    Self Consistent Field (SCF#1)
===================================================================
Iteration     Free Energy (Ha/atom)   SCF Error        Timing (sec)
1            -6.0722904791E+00        3.403E-01        0.631
2            -6.1407156894E+00        7.551E-01        0.273
3            -6.0341313172E+00        7.503E-02        0.258
4            -6.0307338577E+00        1.887E-02        0.245
5            -6.0305973037E+00        9.052E-03        0.253
6            -6.0305318908E+00        2.729E-03        0.246
7            -6.0305402408E+00        2.668E-03        0.252
8            -6.0305483630E+00        2.218E-03        0.244
9            -6.0305647823E+00        1.691E-03        0.239
10           -6.0305675872E+00        2.620E-03        0.185
11           -6.0305690822E+00        2.936E-03        0.223
12           -6.0305900698E+00        3.748E-03        0.222
13           -6.0306038272E+00        3.957E-03        0.232
14           -6.0306148810E+00        2.907E-03        0.229
15           -6.0306312091E+00        2.023E-03        0.216
16           -6.0306455946E+00        2.393E-03        0.194
17           -6.0306469719E+00        2.634E-03        0.209
18           -6.0306458159E+00        2.112E-03        0.205
19           -6.0306486144E+00        1.768E-03        0.211
20           -6.0306582464E+00        1.536E-03        0.203
21           -6.0306578714E+00        2.154E-03        0.205
22           -6.0306566894E+00        1.092E-03        0.197
23           -6.0306580555E+00        4.574E-04        0.197
Total number of SCF: 23
====================================================================
                    Energy and force calculation
====================================================================
Free energy per atom               : -6.0306580555E+00 (Ha/atom)
Total free energy                  : -1.8091974166E+01 (Ha)
Band structure energy              : -3.8128622730E+00 (Ha)
Exchange correlation energy        : -4.8662588158E+00 (Ha)
Self and correction energy         : -2.6465053673E+01 (Ha)
-Entropy*kb*T                      : -1.2438389484E-09 (Ha)
Fermi level                        : -2.9120516509E-01 (Ha)
RMS force                          :  7.5281296705E-02 (Ha/Bohr)
Maximum force                      :  1.1195170156E-01 (Ha/Bohr)
Time for force calculation         :  0.009 (sec)
Pressure                           : -3.5253797930E+00 (GPa)
Maximum stress                     :  5.1281197310E+00 (GPa)
Time for stress calculation        :  0.022 (sec)
===================================================================
                    Self Consistent Field (SCF#2)
===================================================================
Iteration     Free Energy (Ha/atom)   SCF Error        Timing (sec)
1            -6.0403347565E+00        2.197E-01        0.303
2            -6.0145170193E+00        1.350E-01        0.260
3            -6.0106136045E+00        7.992E-02        0.248
4            -6.0095487310E+00        3.264E-02        0.254
5            -6.0091559402E+00        5.406E-03        0.249
6            -6.0092081693E+00        3.435E-03        0.250
7            -6.0092493507E+00        3.072E-03        0.243
8            -6.0092729372E+00        1.679E-03        0.250
9            -6.0092798176E+00        2.235E-03        0.219
10           -6.0092830809E+00        3.400E-03        0.203
11           -6.0092788918E+00        3.860E-03        0.191
        """
        with open(out_file, "w") as fd:
            fd.write(content)
        out_dict = _read_out(out_file)
        print(out_dict)
        steps = out_dict["out"]["ionic_steps"]
        assert len(steps) == 2
        assert "total free energy" not in steps[1]
        assert len(steps[1]["convergence"]) == 0


def test_output_parser_all():
    from sparc.common import repo_dir
    from sparc.sparc_parsers.out import _read_out

    for f_out in test_output_dir.glob("**/*.out"):
        data_dict = _read_out(f_out)
