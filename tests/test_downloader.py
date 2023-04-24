import pytest


def test_download(fs, monkeypatch):
    from sparc import download_data
    from sparc.download_data import download_psp

    fs.create_dir("fake")
    download_psp(psp_dir="fake")
