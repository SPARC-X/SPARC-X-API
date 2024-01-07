import os
from pathlib import Path

import pytest

curdir = Path(__file__).parent
psp_dir = curdir / "psps"


def test_parse_psp8_header_valid():
    """Test valid psp8 header"""
    from sparc.sparc_parsers.pseudopotential import parse_psp8_header

    text = """Na    ONCVPSP-4.0.1  r_core=   1.45370   2.30420
     11.0000      9.0000      200312    zatom,zion,pspd
     8       2   1     4   600     0    pspcod,pspxc,lmax,lloc,mmax,r2well
  5.99000000  0.00000000  0.00000000    rchrg fchrg qchrg
     2     2     0     0     0    nproj
     1     1           extension_switch
   0                         2.8120860445295E+00 -1.1184174774003E+00
     1  0.0000000000000E+00 -1.2720057229743E-09  1.5996190483314E-10
    """
    psp_data = parse_psp8_header(text)
    assert psp_data["symbol"] == "Na"
    assert psp_data["zatom"] == 11.0
    assert psp_data["zion"] == 9.0
    assert psp_data["pspxc"] == 2


def test_parse_psp8_header_invalid():
    """Test invalid psp8 header"""
    from sparc.sparc_parsers.pseudopotential import NotPSP8Format, parse_psp8_header

    text = """Cl    ONCVPSP-4.0.1  r_core=   1.45370   2.30420
     11.0000      9.0000      200312    zatom,zion,pspd
     8       2   1     4   600     0    pspcod,pspxc,lmax,lloc,mmax,r2well
  5.99000000  0.00000000  0.00000000    rchrg fchrg qchrg
     2     2     0     0     0    nproj
     1     1           extension_switch
   0                         2.8120860445295E+00 -1.1184174774003E+00
     1  0.0000000000000E+00 -1.2720057229743E-09  1.5996190483314E-10
    """
    with pytest.raises(NotPSP8Format):
        psp_data = parse_psp8_header(text)


def test_parse_psp8_wrong_format():
    """Test a upf header"""
    from sparc.sparc_parsers.pseudopotential import NotPSP8Format, parse_psp8_header

    text = """<UPF version="2.0.1">
  <PP_INFO>

 Ag_ONCV_PBE_FR-1.0.upf

 This pseudopotential file has been produced using the code
 ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)
 fully-relativistic version 2.1.1, 03/26/2014 by D. R. Hamann
 The code is available through a link at URL www.mat-simresearch.com.
 Documentation with the package provides a full description of the
 input data below.

 While it is not required under the terms of the GNU GPL, it is
 suggested that you cite D. R. Hamann, Phys. Rev. B 88, 085117 (2013)
 in any publication using these pseudopotentials.

 -----
    """
    with pytest.raises(NotPSP8Format):
        psp_data = parse_psp8_header(text)


def test_pseudo_infer():
    """Test pseudopotential inferring"""
    from sparc.sparc_parsers.pseudopotential import (
        MultiplePseudoPotentialFiles,
        NoMatchingPseudopotential,
        infer_pseudo_path,
    )

    with pytest.raises(NoMatchingPseudopotential):
        infer_pseudo_path("As", ".")

    assert infer_pseudo_path("As", psp_dir).name == "33_As_15_1.8_2.1_pbe_n_v1.0.psp8"

    with pytest.raises(NoMatchingPseudopotential):
        infer_pseudo_path("Hf", psp_dir)

    with pytest.raises(MultiplePseudoPotentialFiles):
        infer_pseudo_path("Ag", psp_dir)


def test_pseudo_find_and_mapping():
    """Test psp find with infer or direct mapping"""
    from sparc.sparc_parsers.pseudopotential import (
        MultiplePseudoPotentialFiles,
        NoMatchingPseudopotential,
        find_pseudo_path,
    )

    with pytest.raises(NoMatchingPseudopotential):
        find_pseudo_path("Hg", search_path=None, pseudopotential_mapping={})

    with pytest.raises(NoMatchingPseudopotential):
        find_pseudo_path("Hg", search_path=psp_dir, pseudopotential_mapping={})

    with pytest.raises(MultiplePseudoPotentialFiles):
        find_pseudo_path("Ag", search_path=psp_dir, pseudopotential_mapping={})

    with pytest.raises(NoMatchingPseudopotential):
        find_pseudo_path(
            "Ag", search_path=None, pseudopotential_mapping={"Ag": "Ag-PBE.pot"}
        )

    # Use path instead of name
    assert (
        find_pseudo_path(
            "Ag",
            search_path=None,
            pseudopotential_mapping={"Ag": psp_dir / "Ag-PBE.pot"},
        ).name
        == "Ag-PBE.pot"
    )

    assert (
        find_pseudo_path(
            "Ag",
            search_path=psp_dir,
            pseudopotential_mapping={"Ag": "Ag-PBE.pot"},
        ).name
        == "Ag-PBE.pot"
    )

    assert (
        find_pseudo_path(
            "As",
            search_path=psp_dir,
            pseudopotential_mapping={"Ag": "Ag-PBE.pot"},
        ).name
        == "33_As_15_1.8_2.1_pbe_n_v1.0.psp8"
    )

    # None-existing pseudopotential file should still be ok
    assert (
        find_pseudo_path(
            "Ag",
            search_path=psp_dir,
            pseudopotential_mapping={"Ag": "Ag-PBE-Local.pot"},
        ).name
        == "Ag-PBE-Local.pot"
    )


def test_copy_psp():
    """Copy psp"""
    import tempfile

    from sparc.sparc_parsers.pseudopotential import (
        MultiplePseudoPotentialFiles,
        NoMatchingPseudopotential,
        copy_psp_file,
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        source = psp_dir / "Ag-PBE.pot"
        pot = copy_psp_file(source, tmpdir, use_symbol=False)
        assert pot == "Ag-PBE.pot"

        pot = copy_psp_file(source, tmpdir, use_symbol=True)
        assert pot == "Ag.psp8"

        # Rewrite should be possible
        pot = copy_psp_file(source, tmpdir, use_symbol=True)
        assert pot == "Ag.psp8"
