# Maintaining SPARC C/C++ Code

The documentation for SPARC C/C++ code is already complete as its
own. This maintainers' guide is aimed to cover topics like test cases,
CI/CD management, and package release.

## Installation of SPARC C/C++ code

Please refer to SPARC's
[README](https://github.com/SPARC-X/SPARC?tab=readme-ov-file#2-installation)
for more details.

## Test cases

To run the test cases, please refer to the [main
README](https://github.com/SPARC-X/SPARC?tab=readme-ov-file#4-execution)
and [tests'
README](https://github.com/SPARC-X/SPARC/blob/master/tests/README.md).

The test script `SPARC_testing_script.py` was designed to parse SPARC
outputs without the usage of SPARC-X-API, and could be used as a
reference script for testing SPARC-X-API.

### Socket interface test cases

Currently, test cases for socket communication implemented in C/C++
code are placed under
[`tests/Socket/`](https://github.com/SPARC-X/SPARC/tree/master/tests/Socket),
and are not covered by the `SPARC_testing_script.py` script, as most
of them are dependent on the socket calculator class in ASE.

```{note}
As a rule of thumb, the test cases should not reply on SPARC-X-API's
`sparc.socketio` module (and `sparc.SPARC` calculator when `use_socket=True`),
as the design of the Python API may change in the future.
```

These tests can be run in an environment with ASE and SPARC-X-API (for
file I/O alone) installed. Assume that you have compiled `sparc` with socket support into `<SPARC-root>/lib/sparc`,a minimal setup could look like:
```{code} bash
pip install "ase>=3.23"
# Use minimal SPARC-X-API version
pip install git+https://github.com/SPARC-X/SPARC-X-API.git@v1.0.5
cd tests/Socket
# Setup pseudopotential and SPARC command
export ASE_SPARC_COMMAND="mpirun -n 16 ../lib/sparc"
export SPARC_PSP_PATH="../../psps"
# Run all the tests
bash run_all.sh
```

## CI/CD workflows

Currently there are 2 CI/CD workflows when submit commits / PRs to the
public SPARC C/C++ repo:

1) Building SPARC binary and run unit test [link](https://github.com/SPARC-X/SPARC/blob/master/.github/workflows/build-test.yml)
2) LaTeX documentation and parameters validation [link](https://github.com/SPARC-X/SPARC/blob/master/.github/workflows/update-doc-pdf.yml)

Like all CI/CD workflows [in SPARC-X-API](#cicd-sparc-x-api), the
workflows contain the `workflow_dispatch` keywords, allowing them to
be manually executed / re-run from a personal fork. Please check the
sections below for more details.

### Unit test workflow

This workflow contains two parts
1) Make sure the SPARC version date (defined in `initialization.c`) is
   synchronized with that in the `ChangeLog`. (Checked by [`test-missing-parameters.py`](https://github.com/SPARC-X/SPARC/blob/master/.github/workflows/test-missing-parameters.py))
2) Make sure the SPARC code compiles and runs on simple test cases

The compilation of SPARC C/C++ code uses the `conda-build` recipe under [`.conda/meta.yaml`](https://github.com/SPARC-X/SPARC/blob/master/.conda/meta.yaml), and
Only the `quick_run` test cases are checked during this step
```{code} bash
python SPARC_testing_script.py quick_run
```


### Documentation and parameters validation workflow

This is a quick test to ensure that:

1. The LaTeX documentation under `doc/.LaTeX` and subdirs can be
compiled correctly
2. All the parameter keywords from `.inpt` files
under the `tests` directory have been documented. (Checked by [`test-outdated-package.py`](https://github.com/SPARC-X/SPARC/blob/master/.github/workflows/test-outdated-package.py))

The LaTeX files are compiled using the light-weight [`tectonic`
engine](https://tectonic-typesetting.github.io/en-US/),


## Maintaining the conda-forge release

Compiled binary programs of SPARC C/C++ codes are released in
conda-forge channel under the name
[`sparc-x`](https://anaconda.org/conda-forge/sparc-x).  The release is
managed by
[`sparc-x-feedstock`](https://github.com/conda-forge/sparc-x-feedstock). Please
note that this repository is under the conda-forge organization. If
you wish to become a maintainer, please ping
[@alchem0x2a](https://github.com/alchem0x2A).

The feedstock is set to track [new
releases](https://github.com/SPARC-X/SPARC/releases) in SPARC-X, and in most cases the changes should only be added to [`recipe/meta.yaml`](https://github.com/conda-forge/sparc-x-feedstock/blob/main/recipe/meta.yaml).

If the future release mechanism of SPARC changes, it should be
reflected in both the release tag and the following lines in `recipe/meta.yaml`:
```{code} yaml
source:
  url: https://github.com/SPARC-X/{{ package_name }}/archive/refs/tags/v{{ version }}.zip
  sha256: "full sha256 checksum of the above zip"
```

```{note}
conda-forge does not allow direct download from the github commit. If a minor update (e.g. nightly builds) is to be distributed by conda-forge, please create a release in the SPARC main repo first.
```

Once a new release is created in the SPARC main repo, the auto bot
will send a pull request with the updated recipe (or you can manually
create the PR yourself). After confirming that all the checks have
passed, you can merge the pull request to include the compiled
binaries in the conda-forge channel.

### Debug the recipe using local build

Building the recipe locally is identical to the steps for [`sparc-x-api-feedstock`](#conda-forge-build-locally-api).
You need
both the [docker engine](https://docs.docker.com/engine/) and a [conda distribution](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) installed on your local
machine, and run the following command at the root of the local clone of `sparc-x-feedstock`:
```{code} bash
python build-locally.py
```

### Revised builds

If by any chance you find the newly released SPARC binaries in
conda-forge channel contain errors or missing parts (e.g. not working
properly on one platform, missing a shared library, error from
MPI/Lapack dependencies, etc) due to a bug in `recpie/meta.yaml`,
please update the recipe and only change the build number. For example,
if the current recipe contains following build information:
```{code} yaml
build:
  number: 0
```

After implementing the necessary changes in the recipe, please bump the build number to 1:
```{code} yaml
build:
  number: 1
```

This allows conda-forge to distinguish two builds without affecting
the version number, as the actual package is named like
`<os>-<arch>/sparc-x-2.0.0-<hash>_<build>.conda`.

```{note}
If the error stems from the C/C++ itself, you should update the release in the SPARC github repo instead.
```

### Cross-platform compilation

Settings for the cross-platform compilation are defined in the
[`conda-forge.yml`](https://github.com/conda-forge/sparc-x-feedstock/blob/main/conda-forge.yml),
which can be safely edited by yourself (the feedstock maintainer),
although it is usually not required.

Currently the `sparc-x` package is compiled on linux platform alone
with `x86_64` and `aarch64` variants. To change the settings including
host platforms, compilation toolchains (e.g. use native compilation vs
cross-compilation), please refer to the [conda-forge
manual](https://conda-forge.org/docs/maintainer/conda_forge_yml/).

### Future development
There are
several issues for future releases:
- Allow fetching dated releases of SPARC C/C++ code
- Add more build variants for MKL and MPICH
