import os
from pathlib import Path

import pytest
from packaging import version

curdir = Path(__file__).parent


def test_sparc_api(monkeypatch):
    from sparc.api import SparcAPI
    from sparc.utils import locate_api

    monkeypatch.delenv("SPARC_DOC_PATH", raising=False)
    default_ver = SparcAPI().sparc_version
    # No location provided, use default version
    assert default_ver == locate_api().sparc_version
    # Directly load from another doc.
    # Version not detected since src not presented
    older_ver = locate_api(doc_path=curdir / "sparc-latex-doc-202302").sparc_version
    assert older_ver is None
    # Specify SPARC_DOC_PATH
    monkeypatch.setenv(
        "SPARC_DOC_PATH", (curdir / "sparc-latex-doc-202302").resolve().as_posix()
    )
    older_version = locate_api().sparc_version
    assert older_version is None


def test_sparc_params():
    if "SPARC_DOC_PATH" not in os.environ:
        pytest.skip("No $SPARC_DOC_PATH set. Skip")

    from sparc.utils import locate_api

    # Use the default api with SPARC_DOC_PATH
    api = locate_api()
    if api.sparc_version is None:
        pytest.skip("SPARC version not known. skip")

    if version.parse(api.sparc_version) > version.parse("2023.09.01"):
        assert "NPT_SCALE_VECS" in api.parameters
        assert "NPT_SCALE_CONSTRAINTS" in api.parameters
        assert "TWIST_ANGLE" in api.parameters
