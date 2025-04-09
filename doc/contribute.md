# How to Contribute

## Submitting issues and pull requests
We welcome users of SPARC-X and SPARC-X-API to submit
[issues](https://github.com/SPARC-X/SPARC-X-API/issues) and [pull
requests](https://github.com/SPARC-X/SPARC-X-API/pulls).
When reporting a bug, please make sure to include the following
information:

- `SPARC` version (if available). Should look like "Month Day, Year" in the `.out` file)
- `SPARC-X-API` version or commit hash
- Minimal example for reproducing the error
- Error trace message

## Notes for developers

We recommend the following steps to setup the test environment and modify codes.

(dev-env-setup)=
### Setting up environment

We recommend to work in a clean conda environment (or virtualenv),
with pip installation from [master
branch](https://github.com/SPARC-X/SPARC-X-API) (or your own fork),
and download the pseudopotential files.

```bash
conda create -n sparc-x-api-dev python=3.11 pip
conda activate sparc-x-api-dev
git clone https://github.com/SPARC-X/SPARC-X-API.git
cd SPARC-X-API
pip install -e ".[test]"
python -m sparc.download_data
```

The above steps will also install `pre-commit` as the syntax checker and
cleaner before git commits. Although optional, we strongly recommend enabling
`pre-commit` hooks on your local system to ensure the same formatting across
different contributors:
```bash
# At the SPARC-X-API root
pre-commit install
```
which installs the hook to `.git/hooks/pre-commit`.

Please check
[`.pre-commit-config.yaml`](https://github.com/SPARC-X/SPARC-X-API/blob/master/.pre-commit-config.yaml)
for pre-commit hooks used in this project, and change them if needed.

### Installing SPARC C/C++ binary

If you need to test running DFT using the API, compile or install the
`sparc` executables following the
[manual](https://github.com/SPARC-X/SPARC/blob/master/README.md). Check
[some examples](#install-binary) for our recommended
approaches.


### Running tests

All unit tests are in the `tests/` directory and are based on the
[`pytest`](https://docs.pytest.org/en/stable/) framework.

To run all compatible tests at once:
```python
python -m pytest -svv tests/
```

There are several tests that will only be activated when certain
environment variables are set (may subject to changes):

-  [`test_read_all_examples.py`](https://github.com/SPARC-X/SPARC-X-API/blob/master/tests/test_read_all_examples.py): test for parsing all examples from a SPARC C/C++ release. Requires:
   - `SPARC_TESTS_DIR`: directory to local SPARC test files, e.g. `SPARC-master/tests`
   - `SPARC_DOC_PATH` (optional): directory to LaTeX documentation sources, e.g. `SPARC-master/doc/.LaTeX`, to generate compatible JSON schema
   - Valid SPARC command (e.g. via `ASE_SPARC_COMMAND`): optional, to run `test_quick_examples` step


- [`test_socket.py`](https://github.com/SPARC-X/SPARC-X-API/blob/master/tests/test_socket.py): test for socket calculations. Requires:
  - Setting up correct SPARC command (e.g. via `ASE_SPARC_COMMAND`)
  - Having socket-compatible SPARC binary

### Checking test coverage

You could use the `coverage` package to generate a coverage report for
the test codes. The current coverage report for the master branch of
SPARC-X-API can be accessed [here](test_coverage.md).

If running locally, please use the following commands:
```bash
# Set up proper environment variables first
# Run the code at the repo root
coverage run -a -m pytest -svv tests/
coverage html --omit="tests/*.py"
```

which will generate a folder `htmlcov` under the repo root. Open
`htmlcov/index.html` in a browser to see the coverage broken down to
files and lines, as shown in the following screenshot:
```{figure} img/screenshots/coverage_example.png
:figwidth: 80 %
:align: center
```

(doc-edit)=
### Editing documentation

Source files for documentation are placed under `doc/` directory,
which are written using [MyST
flavor](https://myst-parser.readthedocs.io/en/latest/) of
Markdown. The source `.md` files, together with the main `README.md`
are then rendered to html files using `sphinx`.

After [setting up the test environment](setup_environment.md),
additional packages for doc rendering can be installed via:
```bash
cd SPARC-X-API
pip install -e ".[doc]"
```

To generate the doc files locally, run:
```bash
cd SPARC-X-API
sphinx-build doc doc/_build
```
and then open `doc/_build/index.html` in a browser.

Details about hosting the doc on Github pages please refer to the
[maintenance guide](maintainers.md).

```{note}
1. Do not add contents directly in [`index.md`](https://github.com/SPARC-X/SPARC-X-API/blob/master/doc/index.md) except the `{toctree}` block for editing crosslinks.
2. [`README.md`](https://github.com/SPARC-X/SPARC-X-API/blob/master/README.md) should be kept as concise as possible.
```

### Adding examples

All examples are listed in `examples/` directory. Please add examples
that are important for demonstrating the functionalities of
`SPARC-X-API` while the calculations can be finished using moderate
computating power (e.g. a few minutes with 4 CPU cores).

<!-- The examples can have the name in the format `ex[Number]-[purpose].py`. -->

## Notes for repo maintainers

Please check the [maintenance guide](maintainers.md) for roles
involving admin privileges.

```{toctree}
:maxdepth: 1
test_coverage.md
```
