import os
from packaging import version
from pathlib import Path

import pytest

curdir = Path(__file__).parent


def test_sparc_api():
    from sparc.api import SparcAPI
    from sparc.utils import locate_api

    default_ver = SparcAPI().sparc_version
    # No location provided, use default version
    assert default_ver == locate_api().sparc_version
    # Directly load from another doc.
    # Version not detected since src not presented
    older_ver = locate_api(doc_path=curdir / "sparc-latex-doc-202302").sparc_version
    assert older_ver is None

