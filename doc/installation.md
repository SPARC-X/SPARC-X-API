# Installation

## Installting the API
SPARC-X-API may be installed via the following approaches:

(use-conda)=
### Using [`conda`](https://docs.conda.io/en/latest/) (recommended)

You can use any of [`anaconda`](https://docs.anaconda.com/),
[`miniconda`](https://docs.anaconda.com/miniconda/), or
[`micromamba`](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)
tools to install a `conda` package engine. The installation step varies by the platform and CPU architecture. The following example shows the commands to install the latest `miniconda` on x86_64 Linux:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# Always inspect the installation script content and checksum before installation
more Miniconda3-latest-Linux-x86_64.sh
# The interactive installation involves several yes | no promts
bash Miniconda3-latest-Linux-x86_64.sh
# When finished, start another shell session to activate the base environment
```


The rest of the steps will
be made in a [conda
environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands)
to get the latest release of SPARC-X-API that comes with the
pseudopotentials installed:

```bash
# Change 'sparc-env' to your desired name if needed
conda create -c conda-forge -n sparc-env python=3.11 sparc-x-api
conda activate sparc-env
```


On Linux platforms (x86_64, aarch64), you can also install the
pre-compiled SPARC DFT binaries alongside the API:

```bash
# Note the package name is sparc-x
conda install -c conda-forge sparc-x
# Re-activate to have the env variables effective
conda activate sparc-env
```

```{note}
The SPARC binary code distributed by the official conda-forge channel
is **NOT** compiled with socket support (yet). Please refer to the [manual installation](#install-binary) if you wish to install socket-compatible SPARC binary.
```

```{note}
You may choose [`mamba`](https://github.com/mamba-org/mamba) instead of `conda` as the conda engine for faster dependency resolving and installation.
```

(pypi-install)=
### Install from [PyPI](https://pypi.org/project/sparc-x-api/)

If you prefer using `pip` to install SPARC-X-API, there is also a
[mirror package on PyPI](https://pypi.org/project/sparc-x-api/):
```bash
python -m pip install sparc-x-api
```

The pseudopotential files will also be installed in this approach. If you wish to compile the SPARC C/C++ code, please refer to the [manual installation](#install-binary).

```{note}
Please avoid mixing `conda` and `pip` installations for the same package like `sparc-x-api` in the same environment, as this may lead to unexpected behavior.
```


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

```{note}
Please avoid mixing `conda` and `pip` installations for the same package like `sparc-x-api` in the same environment, as this may lead to unexpected behavior.
```

(install-binary)=
## Manual compilation of the SPARC binary code

To utilize the API for drive SPARC calculations, please following the
[SPARC manual](https://github.com/SPARC-X/SPARC) for compilation and
installation of the SPARC DFT code itself. Manual compilation of the latest SPARC
C/C++ code is relatively straightforward as it depends on only a few standard external libraries.

<!-- We recommend using the [`conda-forge` package](#use-conda) to install -->
<!-- the pre-compiled SPARC binary. -->
The examples shown here
compile the SPARC binary code with the following options:
- `USE_MKL=0`: Use `openblas` instead of `MKL` as numerical library.
- `USE_SCALAPACK=1`: Link SPARC with standalone `scalapack` library.
- `USE_FFTW=1`: Use standalone `fftw` for van der Waals functionals.
- `USE_SOCKET=1`: Enable socket communication layer in SPARC.

### Use `conda` toolchains

In the previously configured `sparc-env` environment, install the
build dependencies of SPARC and compile use the `makefile` provided.

The following process compilers SPARC
with OpenMPI/OpenBLAS/Scalapack toolchains.

```bash
conda activate sparc-env
conda install -c conda-forge \
                 make compilers \
		 fftw=*=mpi_openmpi_* \
		 openblas openmpi scalapack
git clone https://github.com/SPARC-X/SPARC.git
cd SPARC/src
# Always inspect the contents of the makefile before continue
cat makefile
make USE_MKL=0 USE_SCALAPACK=1 USE_FFTW=1 USE_SOCKET=1
cd ../..
```

The compiled binary will be at `SPARC/lib/sparc`, and will run when
`sparc-env` environment is activated.

Run a simple command to make sure the SPARC compilation works:
```bash
mpirun -n 1 SPARC/lib/sparc
```

You should see the usage info printed like following:
```
USAGE:
    mpirun -np <nproc> {SPARCROOT}/lib/sparc -name <filename>

    {SPARCROOT} is the location of the SPARC folder

REQUIRED ARGUMENT:
    -name <filename>
           The filename shared by .inpt file and .ion
           file (without extension)

OPTIONS:
    -h, --help
           Display help (from command line).
    -n <number of Nodes>
    -c <number of CPUs per node>
    -a <number of Accelerators (e.g., GPUs) per node>
    -socket <socket>
            <socket> can be either  <host>:<port> or <unix_socket>:UNIX.
            Note: socket (driver) mode is an experimental feature.

EXAMPLE:

    mpirun -np 8 {SPARCROOT}/lib/sparc -name test

    The example command runs sparc with 8 cores, with input file named
    test.inpt, and ion file named test.ion.

NOTE:
    This is a short description of the usage of SPARC. For a detailed
    discription, refer to the manual online at

        https://github.com/SPARC-X/SPARC/tree/master/doc
```

Note that `-socket` is now an accepted option for the SPARC binary.

### Compiling SPARC on HPC

High Performance Clusters (HPC) machines usually have some specific
requirements about the parallel and numerical library setups. While
conda installation works in most cases, it is often true to compile
SPARC with existing MPI/MKL/BLAS libraries to ensure optimal
performance. The following example shows the compilation with Intel
MKL/MPICH on Georgia Tech's [Pheonix Cluster](https://sites.gatech.edu/ewanparktest/phoenix-cluster/):


```bash
module load git intel-one-api fftw
git clone https://github.com/SPARC-X/SPARC.git
cd SPARC/src
make USE_MKL=1 USE_SCALAPACK=0 USE_FFTW=1 USE_SOCKET=1
```

The compiled binary will be at `SPARC/lib/sparc`, and running it
requires the dependent modules to be loaded at runtime. You can
perform the same quick test in the previous section to confirm the
compilation is successful.


Now head to the [setup tutorial](setup_environment.md) to finetune
your settings.
