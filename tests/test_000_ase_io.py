"""This test file has to be run in the beginning of all other unit tests,
or otherwise the test_import_orderX may fail.

It is recommended that the file name of this test not to be changed, or run this test separatedly.
"""
from pathlib import Path

import pytest

curdir = Path(__file__).parent
test_psp_dir = curdir / "psps"
test_output_dir = curdir / "outputs"


def test_import_order1():
    """import ase before sparc"""
    import ase
    from packaging import version

    if version.parse(ase.__version__) >= version.parse("3.23"):
        pytest.skip("")
    from ase.io.formats import ioformats

    assert "sparc" not in ioformats.keys()
    import sparc

    assert "sparc" in ioformats.keys()
    from ase.io import sparc

    assert hasattr(sparc, "read_sparc")
    assert hasattr(sparc, "write_sparc")


def test_import_order2():
    """import ase after sparc"""
    import ase
    from packaging import version

    if version.parse(ase.__version__) >= version.parse("3.23"):
        pytest.skip("")
    from ase.io.formats import ioformats

    import sparc

    assert "sparc" in ioformats.keys()
    from ase.io import sparc

    assert hasattr(sparc, "read_sparc")
    assert hasattr(sparc, "write_sparc")


def test_sparc_fake_write_exp(monkeypatch):
    """Baseline test. Make a fake write_sparc method
    to makesure the sparc.write_sparc works
    """
    import ase
    from packaging import version

    if version.parse(ase.__version__) >= version.parse("3.23"):
        pytest.skip("")

    def fake_write_sparc(filename, atoms, **kwargs):
        print("I'm the fake writer")
        pass

    from pathlib import Path

    import sparc

    monkeypatch.setattr(sparc, "write_sparc", fake_write_sparc)
    from ase.build import bulk

    from sparc import write_sparc

    al = bulk("Al")
    # Both string and PosixPath should work
    write_sparc("test.sparc", al)
    write_sparc(Path("test.sparc"), al)


def test_sparc_fake_write(monkeypatch):
    """Baseline test. Make a fake write_sparc method
    to makesure the ase.io register works
    """
    import ase
    from packaging import version

    if version.parse(ase.__version__) >= version.parse("3.23"):
        pytest.skip("")

    def fake_write_sparc(filename, atoms, **kwargs):
        print("I'm the fake writer")
        pass

    from pathlib import Path

    from ase.io import sparc as _sparc

    import sparc

    monkeypatch.setattr(sparc, "write_sparc", fake_write_sparc)
    monkeypatch.setattr(_sparc, "write_sparc", fake_write_sparc)
    from ase.build import bulk

    al = bulk("Al")
    al.write("test.sparc")
    al.write(Path("test.sparc"))


def test_sparc_fake_read_exp(monkeypatch, fs):
    """Baseline test. Make a fake read_sparc method
    to makesure the sparc.read_sparc register works
    """
    from pathlib import Path

    import ase
    from packaging import version

    if version.parse(ase.__version__) >= version.parse("3.23"):
        pytest.skip("")

    from ase.io import sparc as _sparc

    import sparc

    def fake_read_sparc(filename, *args, **kwargs):
        print("I'm the fake reader")
        from ase.build import bulk

        return bulk("Al")

    monkeypatch.setattr(sparc, "read_sparc", fake_read_sparc)

    from sparc import read_sparc

    fs.create_dir("test.sparc")
    # With current ase version it is not possible to omit the "sparc" part
    # as it always detects the directory method to bundletrajectory
    atoms = read_sparc("test.sparc", format="sparc")
    assert atoms.get_chemical_formula() == "Al"

    atoms = read_sparc(Path("test.sparc"), format="sparc")
    assert atoms.get_chemical_formula() == "Al"


def test_sparc_fake_read(monkeypatch, fs):
    """Baseline test. Make a fake read_sparc method
    to makesure the ase.io register works
    """
    from pathlib import Path

    import ase
    from packaging import version

    if version.parse(ase.__version__) >= version.parse("3.23"):
        pytest.skip("")

    from ase.io import sparc as _sparc

    import sparc

    def fake_read_sparc(filename, *args, **kwargs):
        print("I'm the fake reader")
        from ase.build import bulk

        return [bulk("Al")]

    monkeypatch.setattr(sparc, "read_sparc", fake_read_sparc)
    monkeypatch.setattr(_sparc, "read_sparc", fake_read_sparc)
    from ase.io import read

    fs.create_dir("test.sparc")
    # With current ase version it is not possible to omit the "sparc" part
    # as it always detects the directory method to bundletrajectory
    atoms = read("test.sparc", index=0, format="sparc")
    print(atoms)
    assert atoms.get_chemical_formula() == "Al"

    atoms = read(Path("test.sparc"), index=0, format="sparc")
    assert atoms.get_chemical_formula() == "Al"


def test_sparc_read_auto(monkeypatch, fs):
    """Same version of the fake read but with automatic format discover"""
    from pathlib import Path

    import ase
    from packaging import version

    if version.parse(ase.__version__) >= version.parse("3.23"):
        pytest.skip("")

    from ase.io import sparc as _sparc

    import sparc

    def fake_read_sparc(filename, *args, **kwargs):
        print("I'm the fake reader")
        from ase.build import bulk

        return [bulk("Al")]

    monkeypatch.setattr(sparc, "read_sparc", fake_read_sparc)
    monkeypatch.setattr(_sparc, "read_sparc", fake_read_sparc)
    from ase.io import read

    fs.create_dir("test.sparc")
    # With current ase version it is not possible to omit the "sparc" part
    # as it always detects the directory method to bundletrajectory
    atoms = read("test.sparc", index=0)
    assert atoms.get_chemical_formula() == "Al"

    atoms = read(Path("test.sparc"), index=0)
    assert atoms.get_chemical_formula() == "Al"


def test_ase_io_filetype(fs):
    """If hacked ase.io.formats.filetype correctly recognized sparc format

    Due to the implementation of ase.io.formats, single file tests should be
    done on non-empty files
    """
    import ase
    from packaging import version

    if version.parse(ase.__version__) >= version.parse("3.23"):
        pytest.skip("")
    from ase.io.formats import filetype

    import sparc

    fs.create_dir("test.sparc")
    fs.create_file("test.sparc/test.ion")
    fs.create_file("test.sparc/test.static")
    fs.create_file("test.sparc/test.geopt")
    fs.create_file("test.sparc/test.aimd")
    with open("test.sparc/test.ion", "w") as fd:
        fd.write("\n")
    with open("test.sparc/test.static", "w") as fd:
        fd.write("\n")
    with open("test.sparc/test.geopt", "w") as fd:
        fd.write("\n")
    with open("test.sparc/test.aimd", "w") as fd:
        fd.write("\n")

    assert filetype("test.sparc") == "sparc"
    assert filetype("test.sparc/test.ion") == "ion"
    assert filetype("test.sparc/test.static") == "static"
    assert filetype("test.sparc/test.geopt") == "geopt"
    assert filetype("test.sparc/test.aimd") == "aimd"
