---
title: 'SPARC-X-API: Versatile Python Interface for Real-space Density Functional Theory Calculations'
tags:
  - Density Functional Theories
  - Atomistic Simulations
  - Python
  - Atomic Simulation Environment
  - Socket Interface
authors:
  - name: Tian Tian
    orcid: 0000-0003-0634-0455
    corresponding: true
    affiliation: "1, 3"
  - name: Lucas R Timmerman
    orcid: 0000-0001-5664-5762
    affiliation: 1
  - name: Shashikant Kumar
    orcid: 0009-0001-5134-1580
    affiliation: 1
  - name: Ben Comer
    orcid: 0000-0002-7528-0049
    affiliation: 1
  - name: Andrew J Medford
    orcid: 0000-0001-8311-9581
    corresponding: true
    affiliation: 1
  - name: Phanish Suryanarayana
    orcid: 0000-0001-5172-0049
    corresponding: true
    affiliation: "1, 2"
affiliations:
  - name: College of Engineering, Georgia Institute of Technology, Atlanta, GA 30332, United States of America
    index: 1
  - name: College of Computing, Georgia Institute of Technology, Atlanta, GA 30332, United States of America
    index: 2
  - name: Department of Chemical and Materials Engineering, University of Alberta, Edmonton AB, T6G 2R3, Canada
    index: 3

date: 05 December 2024
bibliography: paper.bib
---

# Summary

Density Functional Theory (DFT) is the de facto workhorse for
large-scale electronic structure calculations in chemistry and materials science.
While plane-wave DFT implementations remain the most widely used,
real-space DFT provides advantages in handling complex boundary
conditions and scaling to very large systems by allowing for the
efficient use of large-scale supercomputers and linear-scaling methods
that circumvent the cubic scaling bottleneck.  The SPARC-X project
([https://github.com/SPARC-X](https://github.com/SPARC-X)) provides
highly efficient and portable real-space DFT codes
for a wide range of first principle applications, available in both
Matlab (M-SPARC [@xu_m-sparc-1.0_2020; @zhang_m-sparc-2.0_2023]) and C/C++
(SPARC [@xu_sparc-1.0_2021; @zhang_sparc-2.0_2024]). The rapid growth
of SPARC’s feature set has created the need for a fully functional
interface to drive SPARC in high-throughput calculations.  Here we
introduce SPARC-X-API, a Python package designed to bridge the SPARC-X
project with broader computational frameworks. Built on the Atomic
Simulation Environment (ASE [@larsen_ase_2017]) standard, the SPARC-X-API
allows users to handle SPARC file formats and run SPARC calculations
through the same interface as with other ASE-compatible DFT packages.
Beyond standard ASE capabilities, SPARC-X-API provides additional
features including 1) support of SPARC-specific setups, including
complex boundary conditions and unit conversion, 2) a JSON schema
parsed from SPARC's documentation for parameter validation and
compatibility checks, and 3) a comprehensive socket communication layer
derived from the i-PI protocol [@ceriotti_i-pi-1.0_2014;
@kapil_i-pi-2.0_2019] facilitating message passing between low-level C
code and the Python interface.  The goal of the SPARC-X-API is to provide an
easy-to-use interface for users with diverse needs and levels of
expertise, allowing for minimal effort in adapting SPARC to existing
computational workflows, while also supporting developers of advanced
real-space methods.

# Statement of Need

DFT has unarguably become one of the cornerstones of electronic
structure simulations in chemical and materials sciences due to its
simplicity and wide range of applicability.  Among the various
numerical implementations of DFT, the plane-wave pseudopotential
method has gained significant popularity, owing to both its robustness
and the maturity of associated software packages. However, despite
their widespread use, plane-wave methods face several long-standing
challenges, mostly related with the reliance of Fourier transformation
to switch between reciprocal and real-space representations, including
1) establishing efficient schemes on massively parallel computing
environments, 2) overcoming the extensive global communication during
Fourier transformation calculations on very large systems, 3)
developing linear scaling routines [@bowler_order_n_dft_2012] and 4)
handle non-periodic boundary conditions for isolated and semi-finite
systems.  <!-- One --> <!-- long-standing challenge in DFT is to
develop methods that overcome the --> <!-- huge computational cost for
solving the Kohn-Sham equation, which --> <!-- scales cubically with
respect to the system size.  This becomes --> <!-- especially
problematic in massively parallel computing environments, --> <!--
where the extensive global communication required during Fourier -->
<!-- transformations limits the scalability, making it challenging to
--> <!-- efficiently simulate very large systems in plane-wave DFT.
In --> <!-- plane-wave methods, the global nature of the Fourier basis
used limits --> <!-- the ability to achieve linear scaling --> <!--
[@bowler_order_n_dft_2012]. Moreover, the periodic nature of the -->
<!-- Fourier basis enforces the use of periodic boundary conditions,
making --> <!-- the simulation setup of isolated and semi-finite
systems --> <!-- non-straightforward.  --> A compelling alternative to
overcome these limitations is to solve the Kohn-Sham equations using a
finite-difference (FD) approach on real-space grids. The locality of
the FD method makes real-space DFT methods inherently scalable and
paves the way for the development of linearly-scaling solutions to the
Kohn-Sham equations.  Real-space DFT also naturally supports periodic
and Dirichlet boundary conditions, and combinations thereof, allowing
for the flexible treatment of systems in any dimensionality.

In the past few years, the SPARC-X project
([https://github.com/SPARC-X](https://github.com/SPARC-X)) has led
efforts to develop an open-source, real-space DFT code that is both
user-friendly and competitive with state-of-the-art plane-wave
codes. The philosophy of the SPARC-X project is to provide codes that
are highly efficient and portable (i.e., straightforward to install and
use across various computational environments). The codes also seek to
be both user- and developer-friendly to facilitate the
implementation of new algorithms. In line with this, SPARC-X offers
real-space DFT algorithms through two implementations: 1) Matlab-based
M-SPARC [@xu_m-sparc-1.0_2020; @zhang_m-sparc-2.0_2023] for algorithm
prototyping and small-system simulations, with no external
dependencies other than Matlab itself, and 2) C/C++ based SPARC
[@xu_sparc-1.0_2021; @zhang_sparc-2.0_2024] for large-scale production
calculations that can accommodate a wide range of system sizes and
requires only MPI and MKL/BLAS for compilation.  New features of
SPARC include spin-orbit coupling, dispersion
interactions, and advanced exchange-correlation (xc) functionals
[@zhang_sparc-2.0_2024], linear-scaling Spectral Quadrature (SQ)
method [@suryanarayana_sparc_sq_2018], cyclic/helical symmetry
[@sharma_sparc_cyclix_2021], real-space density functional
perturbation theory (DFPT) [@sharma_sparc_dfpt_2023], orbital-free DFT
(ODFT) [@ghosh_sparc_ofdft_2016], and on-the-fly machine-learning force
fields (OTF-MLFF) [@kumar_ofdft_delta_ml_2023;
@timmerman_sparc_mlff_2024; @kumar_sparc_mlff_2024].

The rapid development of SPARC has naturally created a demand for a
fully functional and user-friendly interface that integrates SPARC
smoothly into high-throughput workflows, which also completes the SPARC-X toolkit and complements its philosophy of
usability and portability.
To address this, we introduce the SPARC-X-API, a
Python interface designed to bridge the SPARC code with a wide range
of scientific workflows. The SPARC-X-API builds upon the Python
wrapper originally shipped with SPARC version 1.0
[@xu_sparc-1.0_2021], offering an API compatible with the widely-used
ASE (ASE [@larsen_ase_2017]) standard and updated with the latest
versions of SPARC. With ASE's support for various popular DFT methods,
including both plane-wave (e.g. VASP [@kresse_vasp_1996], Quantum
ESPRESSO [@giannozzi_qe_2017], and Abinit [@gonze_abinit_2020]), and
real-space (e.g. GPAW [@enkovaara_gpaw_1_2011; @mortensen_gpaw_2_2024]
and Octopus [@tancogne_dejean_octopus_2020]) implementations,
SPARC-X-API enables seamless integration of SPARC into existing
workflows, allowing users to incorporate real-space DFT calculations
with minimal adjustments.  The modular design of SPARC-X-API makes it
straightforward to be plugged into complex computational workflows,
for example high-throughput dynamics simulations by i-PI
[@litman_i-pi-3.0_2024] and PLUMED [@bonomi_plumed_2019], as well as
active machine learning frameworks including FineTuna
[@musielewicz_finetuna_2022], powered by state-of-art neural network
interatomic potentials such as FAIR-Chem
([https://github.com/FAIR-Chem/fairchem](https://github.com/FAIR-Chem/fairchem))
and MACE-MP [@ilyes_mace_2023] model series.  A summary of the role
SPARC-X-API in the SPARC-X project is shown in
\autoref{fig:sparc-overview}.  In addition to the capabilities
inherited from ASE, SPARC-X-API seeks to enhance the user experience
in a few key aspects, including 1) supporting SPARC-specific features
in an ASE-compatible API, 2) a parameter validation mechanism based on
SPARC's `LaTeX` documentation, and 3) a versatile socket communication
layer for efficient high-throughput calculations. Details will be
discussed next.

![Overview of SPARC-X-API in the SPARC-X project system
\label{fig:sparc-overview}
](fig/fig_sparc_api_overview.svg){ width=90% }



# Features and Functionalities

The SPARC-X-API is structured as a Python package, `sparc`. A summary of
its key functionalities is provided below; for current detailed
documentation, please refer to the [official
documentation](https://sparc-x.github.io/SPARC-X-API).

## `sparc.io`: File I/O Manipulation

In SPARC and M-SPARC calculations, input information is provided
by two files: a `.inpt` (cell dimensions, boundary conditions,
calculation flags), and a `.ion` file (atomic configurations and
locations to pseudopotential). Depending on the type of calculation,
various output files may be written, such as`.static`, `.geopt` or
`.aimd`. The separation of information across multiple files means
converting ASE `Atoms` objects to SPARC input files or retrieving
energy and forces information from SPARC calculations requires
handling more than just a single file, as is common in most ASE I/O
formats. To manage this, the SPARC-X-API operates on the directory level,
treating each calculation directory as a "SPARC bundle". The
`sparc.io.SparcBundle` class facilitates reading from and writing to
this bundle, ensuring that all necessary input and output files are
properly handled. By default, the SPARC-X-API also copies relevant
pseudopotential files into the calculation directory, making the SPARC
bundle portable across different machines. From version 1.0.7 onwards,
the SPARC-X-API leverages the new features introduced in ASE version 3.23
to register as an external I/O format, allowing reading and writing
SPARC files directly using `ase.io` submodule:

```py
from ase.io import read, write
# 1. Read a SPARC bundle by specifying the `sparc` format
atoms = read("sparc_output_dir", format="sparc")
# 2. Write to a SPARC bundle from aboth object
write("sparc_input_dir", atoms, format="sparc")
```

The SPARC-X-API also supports parsing complex boundary conditions from the
`.inpt` file. The periodic (P) and Dirichlet (D) boundary conditions
are translated into `True` and `False` values, respectively, in the
corresponding `pbc` direction of an `Atoms` object. Standard ASE
objects do not natively support cyclic (C) or helical (H) boundary
conditions that are available in SPARC, so the SPARC-X-API
treats them similarly to Dirichlet boundaries
and stores the original boundary condition information in the `info`
attribute of the atomic object. This ensures that the correct boundary
combinations are preserved when re-writing to SPARC input files.


## `sparc.api`: Parameter Validation

In the ASE ecosystem, default calculator interfaces such as
`FileIOCalculator` do not implement parameter validation, which can
lead to issues such as incorrect parameter settings or incompatibility
when running calculations through ASE. To address this, the SPARC-X-API
introduces a robust parameter validation system using a JSON schema
generated from SPARC’s [LaTeX
documentation](https://github.com/SPARC-X/SPARC/tree/master/doc/.LaTeX). A
JSON schema contains the version of the SPARC software, a list of
input parameters used in `.inpt` and `.ion` files, as well as
supported data types and parameter categories. Validation is handled via the `sparc.api.SparcAPI` class, and includes:

- Verify that the schema is compatible with the version of SPARC binary.
- Convert `.inpt` fields into Python data types.
- Validate input parameters in both string and numerical formats.
- Output help information about specific parameter(s).

Each release of the SPARC-X-API contains a copy of a JSON schema linked
with the latest SPARC release as the default validator, although the
users can select different combination of SPARC versions and
schemas depending on the version they are using. The separation between the
SPARC-X-API and the core SPARC code not only
prevents the need for hard-coding parameter lists into the API, but
also facilitates easier maintenance: the "central truth" of parameters
remains in the SPARC documentation, maintained by the SPARC core
developers, while the SPARC-X-API focuses on providing a user-friendly
interface without being tied to constant updates. This approach maximizes
flexibility and avoids version conflicts between the API and the underlying code.

## `sparc.calculator`: Socket-Communication Calculator Interface

The submodule `sparc.calculator` provides a class `SPARC` as the main
entry point for driving SPARC calculations. This class provides two modes
of operation: 1) a file I/O-based calculator extending the
`ase.calculators.FileIOCalculator` class, and 2) a comprehensive
socket communication layer that allows direct communication between
the Python API and low-level C/C++ code.

In file I/O mode, the SPARC calculator object utilizes the
`sparc.io.SparcBundle` for generating input files and
`sparc.api.SparcAPI` for parameter validation, while the mode of
calculation (single-point, relaxation or molecular dynamics) is
controlled by the input flags. For users transitioning from other DFT
packages and their ASE calculators, the SPARC-X-API is designed to
minimize adaptation effort, but the API is designed to also enable
advanced inputs from expert users. The `SPARC` calculator class
achieves this by supporting two sets
of input parameters: 1) lower-case special parameters that follow
conventions from other ASE DFT calculators (e.g. real-space grid
spacing `h` from GPAW, and exchange-correlation keyword `xc` from
VASP) that use the ASE default Angstrom-eV system, and 2) case-insensitive raw SPARC
input parameters in Bohr-Hartree units for fine-grained control. This
dual approach is designed so that users familiar with other DFT codes
can adopt SPARC with minimal changes to their existing
workflows, while expert users can exert full control.
Basic DFT calculations can be covered by using standard ASE
parameter sets in the SPARC-X-API, as shown by the side-by-side
constructor with VASP and GPAW, using the same
exchange-correlation functional and compatible convergence settings:


```py
#1. Using VASP
from ase.calculators.vasp import Vasp
calc = Vasp(xc="pbe", kpts=(9, 9, 9), ecut=450, ediff=1.e-4)

#2. Using GPAW
from gpaw import GPAW
calc = GPAW(xc="pbe", kpts=(9, 9, 9), h=0.25, convergence={"energy": 1.e-4})

#3. Using SPARC
from sparc.calculator import SPARC
calc = SPARC(xc="pbe", kpts=(9, 9, 9), h=0.25, convergence={"energy": 1.e-4})
```

In high-throughput frameworks requiring thousands of single-point DFT
evaluations, relying on file I/O mode can be inefficient, as
calculations are restarted at each DFT call and the total number of
files may exceed SPARC's default file count limit.  The socket layer in
the SPARC-X-API avoids these limitations by directly communicating with a
long-running SPARC process for updating atomic positions, while
keeping density and orbitals in memory and reducing self-consistent
field (SCF) cycles. While alternative communication methods exist,
such as C-binding approaches seen in GPAW [@mortensen_gpaw_2_2024] and
Psi4 [@smith_psi4_2020], these typically involve complex compilation
and integration steps when installing the Python package. We chose a
socket-based communication layer for its simplicity, which allows for
a clear separation between the Python and SPARC codebases, minimal
modifications to the existing C/C++ code, and ease of installation without
requiring recompilation.

The communication protocol used in the SPARC-X-API socket, referred to as the
SPARC protocol, is based on the i-PI protocol
[@ceriotti_i-pi-1.0_2014; @kapil_i-pi-2.0_2019], which is also adopted
by a wide range of ASE calculators. The SPARC protocol introduces
additional header types and supports binary data transfers via
Python's pickle format.  While SPARC’s C/C++ code maintains compatibility
with the original i-PI standard, the SPARC-X-API leverages the extended
protocol with pickle decoding. The two-tier design offers flexibility
for socket calculations. At its core, the SPARC binary can communicate
directly with any i-PI-compatible server, such as
`ase.calculators.socketio.SocketIOCalculator` in ASE, using the basic
protocol, though this requires careful setup by the user.
However, the SPARC-X-API leverages the
SPARC protocol, which allows the API to internally relay more advanced
 data types to the SPARC
binary, handling object decoding and socket resets automatically. When
running socket calculations on a single machine, users can activate
socket mode by simply adding `use_socket=True` to the `SPARC`
calculator constructor, enabling UNIX socket communication without
additional setup. More importantly, the design of the SPARC protocol
allows easy and seamless integration in distributed computational
systems, offering the following features: 1) flexible client
initialization / restart 2) efficient data transfer 3) heterogeneous
computational setups.
The design of the SPARC protocol allows insertion of bidirectional
additional routines between two DFT calls, allowing further control
over the low-level C/C++ code.
\autoref{fig:socket-hetero} summarizes the
server-client setup across hybrid computing platforms.

![Example of socket communication across hybrid computing platforms using SPARC-X-API
\label{fig:socket-hetero}
](fig/fig_socket_hetero.svg){ width=100% }



## Miscellaneous Helper Functionalities

The SPARC-X-API also provides several helper functions to facilitate user
installation and testing, including:

- `sparc.quicktest`: a utility to verify the installation and
  environment setups for `SPARC-X-API` and `SPARC`.
- `sparc.docparser`: a submodule to convert existing `LaTeX`
  documentation included in SPARC source code into JSON schema.
- `sparc.download_data`: a tool to download the latest ONCV
  pseudopotentials distributed by SPARC.
- `sparc-ase`: an extension to the commandline `ase` tool, adding
  compatibility with SPARC file formats.

# Code Release and Maintenance

The SPARC-X-API is released as source code in github repository
[https://github.com/SPARC-X/SPARC-X-API](https://github.com/SPARC-X/SPARC-X-API),
and as a `conda-forge` package
[`sparc-x-api`](https://anaconda.org/conda-forge/sparc-x-api). When
installed using `conda-forge`, the package is bundled with the
optimized SPMS pseudopotentials [@shojaei_sparc_pseudopot_2023], and
compatible with the
[`sparc`](https://anaconda.org/conda-forge/sparc-x) package that
contains the compiled SPARC binary.

It also integrates continuous integration (CI) workflows for:

- Unit testing and code coverage
- Fetching the latest SPARC documentation for updating the JSON schema
- Validating all test examples from the SPARC repository

These workflows ensure that SPARC-X-API remains up-to-date with
ongoing SPARC developments while separating parameter updates from the
main SPARC maintainers’ efforts.

# Acknowledgements

The authors gratefully acknowledge the support of the U.S. Department
of Energy, Office of Science, under Grant No. DE-SC0019410 and
DE-SC0023445.


# References
