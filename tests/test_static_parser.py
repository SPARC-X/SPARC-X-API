import os
import tempfile
from pathlib import Path

import numpy as np
import pytest
from ase.units import Bohr, Hartree

curdir = Path(__file__).parent
test_output_dir = curdir / "outputs"
repo_dir = curdir.parent


def test_static_parser():
    from sparc.sparc_parsers.static import _read_static

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        static_file = tmpdir / "test.static"

        with open(static_file, "w") as fd:
            fd.write(
                """***************************************************************************
                            Atom positions
***************************************************************************
Fractional coordinates of Fe:
      0.0000000000       0.0000000000       0.0000000000
      0.5200000000       0.5100000000       0.4900000000
Total free energy (Ha): -2.283157353113279E+02
Atomic forces (Ha/Bohr):
  8.0738249305E-01   3.7399117306E-01  -3.5796157735E-01
 -8.0738249305E-01  -3.7399117306E-01   3.5796157735E-01
Stress (GPa):
 -2.1918863425E+04   1.3932450782E+03  -5.1023512490E+01
  1.3932450782E+03  -2.1975897437E+04  -1.2676410947E+02
 -5.1023512490E+01  -1.2676410947E+02  -2.2380745784E+04
            """
            )
        data_dict = _read_static(static_file)
        assert "static" in data_dict
        static_dict = data_dict["static"][0]
        assert "atoms" in static_dict
        assert tuple(static_dict["atoms"]["symbols"]) == ("Fe", "Fe")
        assert np.isclose(static_dict["free energy"], -2.283157353113279e02 * Hartree)
        assert static_dict["forces"].shape == (2, 3)
        assert np.isclose(
            static_dict["forces"][0, 0], 8.0738249305e-01 * Hartree / Bohr
        )
        assert static_dict["stress"].shape == (6,)


def test_static_parser_missing_fields():
    from sparc.sparc_parsers.static import _add_cell_info, _read_static

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        static_file = tmpdir / "test.static"
        with open(static_file, "w") as fd:
            fd.write(
                """***************************************************************************
                            Atom positions
***************************************************************************
Fractional coordinates of Al:
      0.0000000000       0.0000000000       0.0000000000
Total free energy (Ha): -2.212996080865029E+00
Atomic forces (Ha/Bohr):
  0.0000000000E+00   0.0000000000E+00   0.0000000000E+00
            """
            )
        data_dict = _read_static(static_file)
        assert "static" in data_dict
        static_dict = data_dict["static"][0]
        assert "atoms" in static_dict
        assert "stress" not in static_dict
        assert tuple(static_dict["atoms"]["symbols"]) == ("Al",)
        cell = (
            np.array([[0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]])
            * 5.656854249492380
            * Bohr
        )
        new_data_dict = _add_cell_info([static_dict], cell)[0]

        assert "coord" in new_data_dict["atoms"]
        assert new_data_dict["atoms"]["coord"].shape == (1, 3)
        assert np.isclose(new_data_dict["atoms"]["coord"][0, 0], 0)


def test_static_parser_no_atoms():
    from sparc.sparc_parsers.static import _add_cell_info, _read_static

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        static_file = tmpdir / "test.static"
        with open(static_file, "w") as fd:
            fd.write(
                """
Total free energy (Ha): -2.212996080865029E+00
            """
            )
        data_dict = _read_static(static_file)
        assert "static" in data_dict
        static_dict = data_dict["static"][0]
        assert "atoms" not in static_dict
        assert "stress" not in static_dict
        cell = (
            np.array([[0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]])
            * 5.656854249492380
            * Bohr
        )
        new_data_dict = _add_cell_info([data_dict], cell)[0]

        assert "atoms" not in new_data_dict


def test_static_multi_image_same_cell():
    """Read a static file with multiple images (like in socket mode)

    In this test, the lattice for each image is the same (ionic relaxation)
    """
    from sparc.sparc_parsers.static import _add_cell_info, _read_static

    cell = np.array([[4.05, 0.0, 0.0], [0.0, 4.05, 0.0], [0.0, 0.0, 4.05]])

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        static_file = tmpdir / "test.static"
        with open(static_file, "w") as fd:
            fd.write(
                """***************************************************************************
                        Atom positions (socket step 1)
***************************************************************************
Fractional coordinates of Al:
      0.0245290938       0.9931721333       0.0319846198
      0.0752113506       0.4884368716       0.4884376815
      0.5779858173       0.0378980123       0.4768160790
      0.5267930889       0.4771151753       0.9770010000
Lattice (Bohr):
  7.6533908096E+00   0.0000000000E+00   0.0000000000E+00
  0.0000000000E+00   7.6533908096E+00   0.0000000000E+00
  0.0000000000E+00   0.0000000000E+00   7.6533908096E+00
Total free energy (Ha): -9.043005588833887E+00
Atomic forces (Ha/Bohr):
  4.5944553076E-03  -2.9203865816E-04  -9.9557951003E-03
 -2.3316843645E-03   5.5280429722E-03   4.4777314303E-03
 -4.0851154286E-03  -1.3853146416E-02   5.0853960964E-03
  1.8223444855E-03   8.6171421021E-03   3.9266757367E-04
Stress (GPa):
 -2.1676825726E+01  -2.6949376288E-01   8.3753040021E-01
 -2.6949376288E-01   1.0418218149E+01   4.4819979749E-03
  8.3753040021E-01   4.4819979749E-03  -2.5633685934E+01
***************************************************************************
                        Atom positions (socket step 2)
***************************************************************************
Fractional coordinates of Al:
      0.0246648938       0.9931635012       0.0316903531
      0.0751424321       0.4886002642       0.4885700296
      0.5778650741       0.0374885506       0.4769663901
      0.5268469531       0.4773698741       0.9770126049
Lattice (Bohr):
  7.6533908096E+00   0.0000000000E+00   0.0000000000E+00
  0.0000000000E+00   7.6533908096E+00   0.0000000000E+00
  0.0000000000E+00   0.0000000000E+00   7.6533908096E+00
Total free energy (Ha): -9.043161352197208E+00
Atomic forces (Ha/Bohr):
  4.3923550661E-03  -2.4509690104E-04  -9.8419548499E-03
 -2.2709031436E-03   5.4167018956E-03   4.4134407123E-03
 -3.9970707817E-03  -1.3709267553E-02   4.9648744789E-03
  1.8756188591E-03   8.5376625582E-03   4.6363965864E-04
Stress (GPa):
 -2.1692356285E+01  -2.6724550796E-01   8.2319727453E-01
 -2.6724550796E-01   1.0434880328E+01   3.2822223336E-03
  8.2319727453E-01   3.2822223336E-03  -2.5491008099E+01
***************************************************************************
                        Atom positions (socket step 3)
***************************************************************************
Fractional coordinates of Al:
      0.0294363605       0.9928972494       0.0209989185
      0.0726755235       0.4944844938       0.4933644049
      0.5735230074       0.0225960074       0.4823597926
      0.5288844593       0.4866444395       0.9775162642
Lattice (Bohr):
  7.6533908096E+00   0.0000000000E+00   0.0000000000E+00
  0.0000000000E+00   7.6533908096E+00   0.0000000000E+00
  0.0000000000E+00   0.0000000000E+00   7.6533908096E+00
Total free energy (Ha): -9.047001594631228E+00
Atomic forces (Ha/Bohr):
  2.4156119100E-03   1.2164081954E-03  -6.0736428757E-03
 -2.4501250186E-03   2.5307557367E-03   2.0646849306E-03
 -2.2623980616E-03  -8.4907372223E-03   2.2296183039E-03
  2.2969111703E-03   4.7435732902E-03   1.7793396413E-03
Stress (GPa):
 -2.0388752465E+01  -1.3826995208E-01   3.3893431251E-01
 -1.3826995208E-01   1.0147311461E+01  -6.6073109432E-03
  3.3893431251E-01  -6.6073109432E-03  -2.2513420091E+01
                """
            )
        data_dict = _read_static(static_file)
        assert "static" in data_dict
        assert len(data_dict["static"]) == 3
        steps = data_dict["static"]
        for step in steps:
            assert "atoms" in step
            assert "coord_frac" in step["atoms"]
            assert "lattice" in step
            assert "free energy" in step
            assert "forces" in step
            assert "stress" in step

            # The lattices are the same
            assert np.isclose(step["lattice"], cell, 1.0e-6).all()


def test_static_multi_image_diff_cell():
    """Read a static file with multiple images (like in socket mode)

    In this test, the lattice for each image is different (ionic relaxation)
    """
    from sparc.sparc_parsers.static import _add_cell_info, _read_static

    cell = np.array([[4.05, 0.0, 0.0], [0.0, 4.05, 0.0], [0.0, 0.0, 4.05]])

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        static_file = tmpdir / "test.static"
        with open(static_file, "w") as fd:
            fd.write(
                """***************************************************************************
                        Atom positions (socket step 1)
***************************************************************************
Fractional coordinates of Al:
      0.0435568480       0.0098804249       0.0241663700
      0.0553306963       0.5461125430       0.4758696820
      0.5234589733      -0.0037372150       0.4974513864
      0.5101382346       0.5035566314       0.0359079878
Lattice (Bohr):
  7.6533908096E+00   0.0000000000E+00   0.0000000000E+00
  0.0000000000E+00   7.6533908096E+00   0.0000000000E+00
  0.0000000000E+00   0.0000000000E+00   7.6533908096E+00
Total free energy (Ha): -9.041897494954259E+00
Atomic forces (Ha/Bohr):
 -6.3263660275E-03   4.1224438525E-03  -9.7573141360E-03
 -1.1198940824E-02  -9.4391851504E-03   1.3724648200E-02
  7.1551844744E-03   6.8108725330E-03   8.1322757442E-03
  1.0370122377E-02  -1.4941312350E-03  -1.2099609809E-02
Stress (GPa):
 -1.2943266427E+01  -7.3004908329E-01   1.0042725753E+00
 -7.3004908329E-01  -1.3976465331E+01   5.9432478434E-01
  1.0042725753E+00   5.9432478434E-01  -8.2525496305E+00
***************************************************************************
                        Atom positions (socket step 2)
***************************************************************************
Fractional coordinates of Al:
      0.0401072929      -0.0151050966      -0.0130412778
     -0.0264930524       0.5213680896       0.4431718840
      0.5430817720      -0.0187952321       0.5078775085
      0.4938427068       0.5361014305      -0.0508676718
Lattice (Bohr):
  7.7299247177E+00   0.0000000000E+00   0.0000000000E+00
  0.0000000000E+00   7.7299247177E+00   0.0000000000E+00
  0.0000000000E+00   0.0000000000E+00   7.7299247177E+00
Total free energy (Ha): -9.043789471327491E+00
Atomic forces (Ha/Bohr):
 -4.5160329928E-03   1.0836654669E-02  -2.0034847850E-03
  7.8472069094E-03  -1.0356602905E-02   4.6601556005E-03
 -5.6452194926E-03   1.4120642271E-02  -9.1708495469E-03
  2.3140455760E-03  -1.4600694035E-02   6.5141787315E-03
Stress (GPa):
 -7.3542604262E+00   1.8630500028E+00  -7.3044231388E-01
  1.8630500028E+00  -2.2936813817E+01   1.6489555596E+00
 -7.3044231388E-01   1.6489555596E+00  -9.3987769886E-01
***************************************************************************
                        Atom positions (socket step 3)
***************************************************************************
Fractional coordinates of Al:
     -0.0102903172      -0.0013893044      -0.0527455826
      0.0405005138       0.4557176399       0.4792161144
      0.5124168251      -0.0307478540       0.4738777230
      0.4775553675       0.5136161492       0.0565977287
Lattice (Bohr):
  7.8064586258E+00   0.0000000000E+00   0.0000000000E+00
  0.0000000000E+00   7.8064586258E+00   0.0000000000E+00
  0.0000000000E+00   0.0000000000E+00   7.8064586258E+00
Total free energy (Ha): -9.049938894899322E+00
Atomic forces (Ha/Bohr):
  3.6109364340E-03  -2.6690542034E-03  -4.0562261364E-03
 -1.3775701683E-02   4.3698977978E-03   5.4623877892E-03
  2.2823680670E-03   5.1324965010E-03   4.1411554263E-03
  7.8823971824E-03  -6.8333400954E-03  -5.5473170791E-03
Stress (GPa):
  5.3653952236E+00   5.5029845435E-01   1.4612463857E+00
  5.5029845435E-01   8.7131123988E-01  -8.1756988882E-01
  1.4612463857E+00  -8.1756988882E-01  -3.3530577155E+01
            """
            )
        data_dict = _read_static(static_file)
        assert "static" in data_dict
        assert len(data_dict["static"]) == 3
        steps = data_dict["static"]
        ratios = [1.0, 1.01, 1.02]
        for i, step in enumerate(steps):
            assert "atoms" in step
            assert "coord_frac" in step["atoms"]
            assert "lattice" in step
            assert "free energy" in step
            assert "forces" in step
            assert "stress" in step

            # The lattices are the same
            assert np.isclose(step["lattice"], cell * ratios[i], 1.0e-6).all()


def test_static_incomplete():
    """Read a multi-image static file that is broken
    due to abrupt stop

    In this test, the lattice for each image is different (ionic relaxation)
    """
    from sparc.sparc_parsers.static import _add_cell_info, _read_static

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        content = """***************************************************************************
                        Atom positions (socket step 1)
***************************************************************************
Fractional coordinates of Al:
      0.0435568469       0.0098804247       0.0241663704
      0.0553306963       0.5461125432       0.4758696815
      0.5234589728      -0.0037372148       0.4974513852
      0.5101382346       0.5035566321       0.0359079877
Lattice (Bohr):
  7.6533908096E+00   0.0000000000E+00   0.0000000000E+00
  0.0000000000E+00   7.6533908096E+00   0.0000000000E+00
  0.0000000000E+00   0.0000000000E+00   7.6533908096E+00
Total free energy (Ha): -9.041897494881022E+00
Atomic forces (Ha/Bohr):
 -6.3263659031E-03   4.1224440526E-03  -9.7573144850E-03
 -1.1198940888E-02  -9.4391853227E-03   1.3724648407E-02
  7.1551844772E-03   6.8108726360E-03   8.1322760704E-03
  1.0370122314E-02  -1.4941313660E-03  -1.2099609992E-02
Stress (GPa):
 -1.2943266261E+01  -7.3004908368E-01   1.0042725771E+00
 -7.3004908368E-01  -1.3976465176E+01   5.9432479567E-01
  1.0042725771E+00   5.9432479567E-01  -8.2525500281E+00
***************************************************************************
                        Atom positions (socket step 2)
***************************************************************************
Fractional coordinates of Al:
      0.0401072938      -0.0151050963      -0.0130412790
     -0.0264930519       0.5213680889       0.4431718840
      0.5430817728      -0.0187952321       0.5078775086
      0.4938427062       0.5361014296      -0.0508676716
Lattice (Bohr):
  7.6533908096E+00   0.0000000000E+00   0.0000000000E+00
  0.0000000000E+00   7.6533908096E+00   0.0000000000E+00
  0.0000000000E+00   0.0000000000E+00   7.6533908096E+00
Total free energy (Ha): -9.039359835390405E+00
Atomic forces (Ha/Bohr):
 -5.6420435958E-03   1.2313126326E-02  -2.8998717685E-03
  9.5779301732E-03  -1.1726647745E-02   5.8218036947E-03
 -7.0526004868E-03   1.5921660989E-02  -1.0900664252E-02
  3.1167139094E-03  -1.6508139571E-02   7.9787323255E-03
Stress (GPa):
 -1.1079730116E+01   2.1369400165E+00  -5.5506801999E-01
  2.1369400165E+00  -2.4342715311E+01   1.9045396297E+00
 -5.5506801999E-01   1.9045396297E+00  -3.7041195294E+00
***************************************************************************
                        Atom positions (socket step 3)
***************************************************************************
Fractional coordinates of Al:
     -0.0102903160      -0.0013893037      -0.0527455827
      0.0405005136       0.4557176395       0.4792161136
      0.5124168247      -0.0307478543       0.4738777235
      0.4775553679       0.5136161481       0.0565977284
Lattice (Bohr):
  7.6533908096E+00   0.0000000000E+00   0.0000000000E+00
  0.0000000000E+00   7.6533908096E+00   0.0000000000E+00
  0.0000000000E+00   0.0000000000E+00   7.6533908096E+00
        """
        static_file = tmpdir / "test.static"
        with open(static_file, "w") as fd:
            fd.write(content)
        data_dict = _read_static(static_file)
        steps = data_dict["static"]
        assert len(steps) == 3
        assert "free energy" not in steps[2]
        assert "forces" not in steps[2]
        assert "stress" not in steps[2]
