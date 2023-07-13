## Notes for developers
### Setting up development environment

We recommend use the pip installation from github
```python
git clone https://github.com/SPARC-X/sparc-dft-api.git
pip install -e "sparc-dft-api[test]"
python -m sparc.download_data
```

which will setup sparc-dft-api along with `pytest` toolchains. 

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
finished using moderate computating power (e.g. a few minutes with 4 cores).