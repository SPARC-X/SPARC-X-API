"""Unit test for order of importing.
Submodules like `docparser` `api` and `download_data` should be independent of ase / numpy,
should such modules are not yet available during processes like conda-forge build
"""
import pytest
import sys

def test_download_data():
    import ase
    original_ase = sys.modules.get('ase')
    sys.modules['ase'] = None
    with pytest.raises(ImportError):
        import ase
    from sparc.download_data import download_psp
    # Recover sys ase
    sys.modules['ase'] = original_ase

def test_api():
    import ase
    original_ase = sys.modules.get('ase')
    sys.modules['ase'] = None
    with pytest.raises(ImportError):
        import ase
    from sparc.api import SparcAPI
    # Recover sys ase
    sys.modules['ase'] = original_ase

def test_docparser():
    import ase
    original_ase = sys.modules.get('ase')
    sys.modules['ase'] = None
    with pytest.raises(ImportError):
        import ase
    from sparc.docparser import SPARCDocParser
    # Recover sys ase
    sys.modules['ase'] = original_ase

def test_normal():
    import ase
    original_ase = sys.modules.get('ase')
    sys.modules['ase'] = None
    with pytest.raises(ImportError):
        import ase
    with pytest.raises(ImportError):
        from sparc.io import read_sparc
    # Recover sys ase
    sys.modules['ase'] = original_ase
