import pytest
import os
from pathlib import Path


def test_run_quicktest():
    """Just import and run quicktest"""
    from sparc.quicktest import main

    main()
