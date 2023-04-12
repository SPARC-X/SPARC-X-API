import pytest

def test_import_order1():
    """import ase before sparc
    """
    import ase
    from ase.io.formats import ioformats
    assert "sparc" not in ioformats.keys()
    import sparc
    assert "sparc" in ioformats.keys()
    from ase.io import sparc
    assert hasattr(sparc, "read_sparc")
    assert hasattr(sparc, "write_sparc")


def test_import_order2():
    """import ase after sparc
    """
    import sparc
    import ase
    from ase.io.formats import ioformats
    assert "sparc" in ioformats.keys()
    from ase.io import sparc
    assert hasattr(sparc, "read_sparc")
    assert hasattr(sparc, "write_sparc")

def test_sparc_fake_write(monkeypatch):
    """Baseline test. Make a fake write_sparc method
    to makesure the ase.io register works
    """
    def fake_write_sparc(atoms, filename, **kwargs):
        pass
    import sparc
    monkeypatch.setattr(sparc, "write_sparc", fake_write_sparc)
    from ase.build import bulk
    al = bulk("Al")
    al.write("test.sparc")


def test_sparc_fake_read(monkeypatch, fs):
    """Baseline test. Make a fake read_sparc method
    to makesure the ase.io register works
    """
    import sparc
    from ase.io import read
    def fake_read_sparc(filename, **kwargs):
        from ase.build import bulk
        return bulk("Al")
    fs.create_dir("test.sparc")
    atoms = read("test.sparc")
    assert atoms.get_chemical_formula() == "Al"
    
