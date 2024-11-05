# Installation

The Python API may be installed via either of the following approaches:

### 1. Via `anaconda` or `miniconda` (recommended)

Set up a [conda environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands) and install the Python API,
which includes the pseudopotential files:

```bash
# Change 'sparc-env' to your desired name if needed
conda create -n sparc-env
conda activate sparc-env
conda install -c conda-forge sparc-x-api
```

On Linux platforms (x86_64, aarch64), you can also install the
precompiled `sparc` DFT binaries alongside the API:

```bash
conda install -c conda-forge sparc-x
conda activate sparc-env   # Re-activate to have the env variables effective
```

### 2. Manual installation from source with `pip`


```bash
python -m pip install git+https://github.com/SPARC-X/SPARC-X-API
```

Optionally, you can download the latest SPMS pseudopotentials and unpacks the pseudopotential files into `<python-lib-root>/site-packages/sparc/psp`:

```bash
python -m sparc.download_data
```


To utilize the API for drive SPARC calculations, please
following the [SPARC manual](https://github.com/SPARC-X/SPARC) for
compilation and installation of the SPARC DFT code itself.

### Post-installation check

We recommend the users to run a simple test after installation and setup:

```bash
python -m sparc.quicktest
```

A proper setup will display the following sections at the output's conclusion:

<img width="500" alt="image" src="https://github.com/alchem0x2A/SPARC-X-API/assets/6829706/95cb712e-4c77-4b14-8130-4961e3c50278">

For using the API to parse SPARC input and output files, it's
essential that the "Import" and "JSON API" tests are successful. For
run SPARC calculations, all tests must pass.

Please refer to the [Setting Up the
Environment](#setting-up-the-environment) or guidance on correctly
configuring the environment variables. If you run into further problems, consult our
[Trouble Shooting](doc/troubleshooting.md).
