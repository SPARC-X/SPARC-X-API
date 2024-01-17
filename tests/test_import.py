"""Unit test for order of importing.
Submodules like `docparser` `api` and `download_data` should be independent of ase / numpy,
should such modules are not yet available during processes like conda-forge build
"""
import sys

import ase
import pytest


def test_download_data(monkeypatch):
    monkeypatch.setitem(sys.modules, "ase", None)
    with pytest.raises(ImportError):
        import ase
    from sparc.download_data import download_psp


def test_api(monkeypatch):
    monkeypatch.setitem(sys.modules, "ase", None)
    with pytest.raises(ImportError):
        import ase
    from sparc.api import SparcAPI


def test_docparser(monkeypatch):
    monkeypatch.setitem(sys.modules, "ase", None)
    with pytest.raises(ImportError):
        import ase
    from sparc.docparser import SparcDocParser
