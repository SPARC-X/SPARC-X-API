## Submitting issues and pull requests
We welcome users of SPARC and its python-API to submit issues and pull requests via github.
When reporting a bug, please make sure to include the following information:

- `SPARC` version (if available. Should look like "Month Day, Year" in the `.out` file)
- `sparc-python-api` version or commit hash
- Minimal example for reproducing the error
- Error trace message

## Notes for developers

We recommend the following steps to setup the test environment and modify codes

### Setting up environment

Pip installation from github's master branch (or your own fork), and download
a copy of the latest pseudopotential files.

```python
git clone https://github.com/SPARC-X/sparc-dft-api.git
pip install -e "sparc-dft-api[test]"
python -m sparc.download_data
```

If you need to test running DFT using the API, compile or install the `sparc` executables following the [manual](https://github.com/SPARC-X/SPARC/blob/master/README.md). 


### Running tests

All unit tests are based on `pytest` and inside `tests/` directory. 
To run all tests (no heavy DFT calculations):
```python
python -m pytest -svv tests/
```

If you are on a HPC environment, you can opt to run a comprehensive test suite with DFT calculations:
```python
python -m pytest -svv tests/test_all_dft.py
```

(*Draft, to be implemented later*)

### Adding examples

All examples are listed in `examples/` directory. Please add examples that are important
for demonstrating the functionalities of `sparc-dft-api` while the calculations can be 
finished using moderate computating power (e.g. a few minutes with 4 CPU cores).

The examples can have the name in the format `ex[Number]-[purpose].py`.

### Code structure

Below is a brief overview of the modules in `sparc-dft-api` with simple explanations
```
sparc
├── __init__.py
├── api.py                 # Includes SparcAPI class for parameter validation
├── calculator.py          # Interface to the SPARC DFT code
├── cli.py                 # `sparc-ase` interface
├── common.py              # Definition of common directories
├── docparser.py           # Function and cli interface for parsing the SPARC DFT document
├── download_data.py       # Cli tool to download pseudopotential files
├── io.py                  # Provides `SparcBundle` class, `read_sparc` and `write_sparc` functions
├── quicktest.py           # Cli tool for post-installation sanity check
├── utils.py               # Common utilities
├── psp/                   # Place-holder directory for pseudopotentials (used for `download_data.py`)
├── sparc_json_api         # Directory for maintaining the JSON API
│   └── parameters.json
├── sparc_parsers          # Parsers for individual SPARC in-/output formats
│   ├── __init__.py        
│   ├── aimd.py
│   ├── atoms.py
│   ├── geopt.py
│   ├── inpt.py
│   ├── ion.py
│   ├── out.py
│   ├── pseudopotential.py
│   ├── static.py
│   └── utils.py
```