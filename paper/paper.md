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
    orcid:
    affiliation: [1,2]
  - name: Lucas Timmerman
    orcid:
    affiliation: 1
  - name: Shashikant Kumar
    affiliation: 2
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
    affiliation: [2, 3]
affiliations:
  - name: School of Civil and Environmental Engineering, Georgia Institute of Technology
	index: 1
  - name: School of Chemical and Biomolecular Engineering, Georgia Institute of Technology
    index: 2
  - name: School of Computational Science and Engineering, Georgia Institute of Technology
    index: 3

date: 24 September 2024
bibliography: paper.bib
---

# Summary

Density Functional Theory (DFT) is the de facto gold standard for
electronic structure calculations in chemistry and materials
science. While plane-wave DFT remains the most widely used, real-space
DFT provides advantages in handling complex boundary conditions and
scaling to very large systems. The SPARC-X project
(https://github.com/SPARC-X) has pioneered highly efficient real-space
DFT codes available in both Matlab [@xu_m-sparc-1.0_2020;
@zhang_m-sparc-2.0_2023] and C [@xu_sparc-1.0_2021;
@zhang_sparc-2.0_2024]. However, the specific input formats for SPARC
have often made it challenging for users accustomed to plane-wave DFT
to transition to real-space methods. To address this, we introduce
SPARC-X-API, a Python interface designed to bridge the SPARC-X project
with broader computational frameworks. Built on the atomic simulation
environment (ASE [@larsen_ase_2017]) standard, SPARC-X-API allows
users to handle SPARC file formats, and run SPARC calculations using
the same interface as with other computational packages. SPARC-X-API
provides additional features beyond the standard ASE package,
including 1) support of complex boundary conditions, 2) a JSON schema
for validating and converting calculation parameters, and 3) a
comprehensive calculator interface with advanced socket-communication
support. SPARC-X-API provides a smooth transition for users from
plane-wave DFT, making the access to real-space DFT calculations more
available and flexible for a wider range of users and computational
workflows.

# Statement of Need

Kohn-Sham Density Functional Theory (DFT) has unargubaly become the
cornerstone of electronic simulations in chemical and materials
sciences due to its simplicity and applications across a wide range of
systems. The popularity of DFT over other first-principle methods in
materials simulation largely stems from the simplicity of the
plane-wave pseudopotential implementation, where convergence is
controlled by simply the plane-wave cutoff energy, and solving the
Kohn-Sham equations can be benefited from highly-optimized Fast
Fourier Transform (FFT) packages. While many non-theoretical
researchers may associate DFT exclusively with plane-wave
implementations, this approach has notable limitations. The periodic
nature of the Fourier basis enforces the use of periodic boundary
conditions, making the simulation setup of isolated and semi-finite
systems non-straightforward. Additionally, the global nature of the
Fourier basis causes plane-wave codes to scale poorly with increasing
numbers of parallel processes. A compelling alternative to overcome
these limitations is to solve the Kohn-Sham equations using a
finite-difference approach on real-space grids. Real-space DFT
naturally supports both periodic and Dirichlet boundary conditions,
allowing for the flexible treatment of systems in any
dimensionality. Furthermore, the locality of the finite-difference
grids makes real-space DFT methods inherently scalable, paving the way
for the development of linearly-scaling solutions to the Kohn-Sham
equations.

<!-- Need review @TT 2024.09.30 --> Despite the advantages of
real-space DFT, plane-wave implementations remain dominant in the
field of computational chemistry and materials science, largely due to
the greater accessibility of plane-wave DFT codes and their more
established programmable interfaces. While real-space DFT offers
significant benefits, there are currently few widely used packages
that provide comprehensive real-space DFT capabilities.

The only notable exception has been GPAW
[@mortensen_gpaw_original_2005; @enkovaara_gpaw_1_2011;
@mortensen_gpaw_2_2024], which originally focused on real-space
finite-difference methods. However, in recent years, the development
of GPAW has shifted its focus toward plane-wave implementations
[@@mortensen_gpaw_2_2024], leaving its finite-difference capabilities
underdeveloped and missing key functionality. In contrast, the SPARC-X
project (https://github.com/SPARC-X) has pioneered efforts to develop
an open-source, real-space DFT code that is both user-friendly and
competitive with state-of-the-art plane-wave codes.

SPARC-X offers real-space DFT algorithms through two implementations:
M-SPARC [@xu_m-sparc-1.0_2020; @zhang_m-sparc-2.0_2023] for
prototyping and small-system simulations, and SPARC
[@xu_sparc-1.0_2021; @zhang_sparc-2.0_2024] for large-scale production
calculations that can accommodate a wide range of system
sizes. Although SPARC has demonstrated its computational efficiency
and features a rich set of algorithms, its adoption has been limited
by the lack of a user-friendly interface that can connect the code to
a broader audience of users and computational tools.

To address this, we introduce SPARC-X-API, a Python-based interface
designed to bridge the SPARC-X code with a broader range of scientific
workflows. Built on the Atomic Simulation Environment (ASE
[@larsen_ase_2017]) standard, SPARC-X-API provides seamless file
read/write support for SPARC files and a feature-complete calculator
interface to the SPARC code. With SPARC-X-API, researchers can easily
incorporate real-space DFT into their workflows using familiar tools
and interfaces, making real-space DFT more accessible to a wider range
of users.

<!-- statement of SPARC-X-API v0.1 -->

<!-- SPARC-X-API philosophy -->


# Features and Functionalities

SPARC-X-API offers two key functionalities:

- File I/O: Through the sparc.io submodule, SPARC-X-API implements
  file read/write support for SPARC file formats, including .inpt and
  .ion files. SPARC-X-API operates on the directory level, treating
  each calculation directory as a "SPARC bundle." From version 1.0
  onwards, SPARC-X-API is fully integrated with ASE (version 3.23),
  automatically registering SPARC as an external I/O format.
- Calculator Interface: The sparc.calculator submodule provides a full
  ASE-compatible calculator interface for running SPARC calculations,
  enabling integration with ASE workflows. <!-- IO and socker -->

Unique Features of SPARC-X-API:

1) Support for Bundled File Formats: Unlike typical single-file DFT implementations, SPARC requires both .inpt and .ion files. SPARC-X-API's design simplifies this by reading and writing at the directory level, streamlining the handling of SPARC bundles.

2) JSON Schema for Parameter Validation: SPARC-X-API ensures parameter
consistency through a JSON schema derived from SPARC's LaTeX
documentation. This guarantees compatibility with SPARC's source code,
offering a robust mechanism for validating and converting parameters.
3) Unit Conversions: SPARC-X-API manages the conversion between atomic
units (Hartree, Bohr) used in SPARC and the eV/Å units in ASE.

## Socket-Communication Calculator Interface
SPARC-X-API’s socket-communication layer allows for efficient and
flexible workflows by reducing the overhead of file I/O. This feature
is particularly useful for iterative calculations, such as structural
optimizations and saddle point searches, where traditional file-based
communication can become a bottleneck.

Key advantages:

Efficiency: Eliminates intermediate file I/O by streaming data
directly between processes.  Speed: Enhances performance in iterative
calculations, critical for large-scale simulations.  Flexibility:
Enables real-time modification of calculation parameters without
restarting processes.  SPARC-X-API implements a backward-compatible
i-PI protocol, allowing both low-level and high-level interfacing with
SPARC's DFT code.

# Code Release and Maintenance
SPARC-X-API maximizes accessibility for users by providing streamlined
installation via the conda-forge channel, where the sparc-x-api
package can be installed with default ONCV pseudopotentials. It also
integrates continuous integration (CI) and continuous deployment (CD)
workflows for:

- Unit testing and code coverage
- Fetching the latest SPARC documentation for updating the JSON schema
- Validating all test examples from the SPARC repository

These workflows ensure that SPARC-X-API remains up-to-date with
ongoing SPARC developments while separating parameter updates from the
main SPARC maintainers’ efforts.

# Acknowledgements


# References
