import os
from pathlib import Path

import pytest
from packaging import version

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
    # Specify SPARC_DOC_PATH
    os.environ["SPARC_DOC_PATH"] = (
        (curdir / "sparc-latex-doc-202302").resolve().as_posix()
    )
    older_version = locate_api().sparc_version
    assert older_version is None
