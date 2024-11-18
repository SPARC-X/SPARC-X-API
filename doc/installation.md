# Installation

## Installting the API
SPARC-X-API may be installed via the following approaches:

(use-conda)=
### Using [`conda`](https://docs.conda.io/en/latest/) (recommended)

You can use any of [`anaconda`](https://docs.anaconda.com/),
[`miniconda`](https://docs.anaconda.com/miniconda/), or
[`micromamba`](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)
tools to install a `conda` package engine.  The rest of the steps will
be made in a [conda
environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands)
to get the latest release of SPARC-X-API that comes with the
pseudopotentials installed:

```bash
# Change 'sparc-env' to your desired name if needed
conda create -n sparc-env python=3.11 pip mamba
conda activate sparc-env
conda install -c conda-forge sparc-x-api
```


On Linux platforms (x86_64, aarch64), you can also install the
pre-compiled SPARC DFT binaries alongside the API:

```bash
# Note the package name is sparc-x
conda install -c conda-forge sparc-x
# Re-activate to have the env variables effective
conda activate sparc-env
```

(pypi-install)=
### Install from [PyPI](https://pypi.org/project/sparc-x-api/)

If you prefer using `pip` to install SPARC-X-API, there is also a
[mirror package on PyPI](https://pypi.org/project/sparc-x-api/):
```bash
python -m pip install sparc-x-api
```

The pseudopotential files will also be installed in this approach. If you wish to compile the SPARC C/C++ code, please refer to the [manual installation](#install-binary).


(pip-install)=
### Installing from latest source code

The latest version of the SPARC-X-API can also be installed using `pip` from the source code on Github.

```bash
python -m pip install git+https://github.com/SPARC-X/SPARC-X-API
```


```{note}
The pseudopotential files should be manually downloaded in this case.
```

You can download the latest SPMS pseudopotentials and unpacks the pseudopotential files into `<python-lib-root>/site-packages/sparc/psp`:

```bash
python -m sparc.download_data
```


For developers, please check the [how to
contribute](#setting-up-environment) page for setting up a dev-environment for SPARC-X-API.

(install-binary)=
## Install the SPARC binary code

To utilize the API for drive SPARC calculations, please following the
[SPARC manual](https://github.com/SPARC-X/SPARC) for compilation and
installation of the SPARC DFT code itself.

We recommend using the [`conda-forge` package](#use-conda) to install
the pre-compiled SPARC binary. If you want to compile the latest SPARC
C/C++, it is also straightforward:

### Use `conda` toolchains

In the previously configured `sparc-env` environment, install the
build dependencies and compile.

The following process compilers SPARC
with OpenMPI/OpenBLAS/Scalapack toolchains.

```bash
conda activate sparc-env
mamba install -c conda-forge \
                 make compilers \
				 fftw=*=mpi_openmpi_* \
				 openblas openmpi scalapack
git clone https://github.com/SPARC-X/SPARC.git
cd SPARC/src
make USE_MKL=0 USE_SCALAPACK=1 USE_FFTW=1
```

The compiled binary will be at `SPARC/lib/sparc`, and will run when
`sparc-env` environment is activated.

**TODO** MKL available?

### Compiling SPARC on HPC

High Performance Clusters (HPC) machines usually have some specific
requirements about the parallel and numerical library setups. While
conda installation works in most cases, it is often true to compile
SPARC with existing MPI/MKL/BLAS libraries to ensure optimal
performance. The following example shows the compilation with Intel
MKL/MPICH on Georgia Tech's [Pheonix Cluster](https://sites.gatech.edu/ewanparktest/phoenix-cluster/):

**TODO** make sure modules are correct.
```bash
module load git intel-one-api fftw
git clone https://github.com/SPARC-X/SPARC.git
cd SPARC/src
make USE_MKL=1 USE_SCALAPACK=0 USE_FFTW=1
```

**TODO** add module in script.

The compiled binary will be at `SPARC/lib/sparc`, and running it
requires the dependent modules to be loaded at runtime.


Now head to the [setup tutorial](setup_environment.md) to finetune
your settings.
