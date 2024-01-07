import os
from pathlib import Path

import pytest


def test_run_quicktest():
    """Just import and run quicktest"""
    from sparc.quicktest import main

    main()
