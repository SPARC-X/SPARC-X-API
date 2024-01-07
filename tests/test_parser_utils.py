"""Test file for sparc.sparc_parsers.utils
"""
import pytest


def test_strip_comments():
    from sparc.sparc_parsers.utils import strip_comments

    rawtext = "# Test comment\n   \n# Test 2\nATOM_TYPE:"
    # expected_stripped = ["hello world", "foo bar"]
    # expected_comments = ["this is a comment", "another comment"]
    stripped, comments = strip_comments(rawtext)
    assert len(stripped) == 1
    assert len(comments) == 2  # No white line
    assert stripped[0] == "ATOM_TYPE:"
    assert comments[0] == "Test comment"  # No white spaces


def test_strip_comments2():
    from sparc.sparc_parsers.utils import strip_comments

    lines = """ATOM_TYPE: Cu                 # atom type
ATOMIC_MASS: 63.546           # atomic mass (amu)
PSEUDO_POT: ../../../psps/29_Cu_19_1.7_1.9_pbe_n_v1.0.psp8  # pseudopotential file
N_TYPE_ATOM: 4                # number of atoms of this type
# COORD:                      # Cartesian coordinates (au)
COORD_FRAC:                   # fractional coordinates (in lattice vector basis)
   0.0 0.0 0.0
   0.0 0.5 0.5
   0.5 0.0 0.5
   0.5 0.5 0.0
    """
    stripped, comments = strip_comments(lines)
    assert all(["#" not in line for line in stripped])


def test_bisect_and_strip():
    from sparc.sparc_parsers.utils import bisect_and_strip

    text = "    key : value   \n"
    delimiter = ":"
    expected_data = "key"
    expected_value = "value"
    data, value = bisect_and_strip(text, delimiter)
    assert data == expected_data
    assert value == expected_value


def test_read_block_input():
    from sparc.sparc_parsers.utils import read_block_input, strip_comments

    # Define the input block for testing
    block_lines = """ATOM_TYPE: Cu                 # atom type
ATOMIC_MASS: 63.546           # atomic mass (amu)
PSEUDO_POT: ../../../psps/29_Cu_19_1.7_1.9_pbe_n_v1.0.psp8  # pseudopotential file
N_TYPE_ATOM: 4                # number of atoms of this type
# COORD:                      # Cartesian coordinates (au)
COORD_FRAC:                   # fractional coordinates (in lattice vector basis)
   0.0 0.0 0.0
   0.0 0.5 0.5
   0.5 0.0 0.5
   0.5 0.5 0.0
    """
    block, comments = strip_comments(block_lines)

    # Run 1: without validator
    output = read_block_input(block, validator=None)
    assert all(
        [
            key in output.keys()
            for key in [
                "ATOM_TYPE",
                "ATOMIC_MASS",
                "PSEUDO_POT",
                "N_TYPE_ATOM",
                "COORD_FRAC",
            ]
        ]
    )

    # Run 2: Call the function with the input block and mock validator
    class MockValidator:
        parameters = {"ATOM_TYPE": str, "N_TYPE_ATOM": int}

        def convert_string_to_value(self, key, value):
            # Mocking conversion to int
            try:
                value = int(value)
            except Exception:
                pass
            return value

    with pytest.warns(UserWarning):
        output2 = read_block_input(block, validator=MockValidator())
    assert all(
        [
            key in output2.keys()
            for key in [
                "ATOM_TYPE",
                "ATOMIC_MASS",
                "PSEUDO_POT",
                "N_TYPE_ATOM",
                "COORD_FRAC",
            ]
        ]
    )
    assert output2["N_TYPE_ATOM"] == 4


def test_make_reverse_mapping():
    import numpy as np

    from sparc.sparc_parsers.utils import make_reverse_mapping

    # By applying the reverse mapping function twice, it should go back original
    for i in range(5):
        size = np.random.randint(1, 20)
        sort = np.arange(size)
        np.random.shuffle(sort)
        reverse = make_reverse_mapping(sort)
        rere = make_reverse_mapping(reverse)
        assert np.isclose(sort, rere).all(), f"Reverse error! {sort}, {reverse}, {rere}"
