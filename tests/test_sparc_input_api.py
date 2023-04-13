import pytest


def test_sparc_api():
    from sparc.inputs import SparcInputs
    from sparc.inputs import default_json_api

    # the default api should always exist
    assert default_json_api.is_file()
    sis = SparcInputs()
    assert hasattr(sis, "sparc_version")
    assert hasattr(sis, "categories")
    assert hasattr(sis, "parameters")
    assert isinstance(sis.parameters, dict)
    assert hasattr(sis, "other_parameters")


def test_api_validate():
    from sparc.inputs import SparcInputs
    import numpy as np

    sis = SparcInputs()
    # Integer
    assert sis.validate_input("CALC_PRES", "0")
    assert sis.validate_input("CALC_PRES", "0 ")
    assert sis.validate_input("CALC_PRES", 0)
    assert sis.validate_input("CALC_PRES", True)
    # Invalid parameters
    with pytest.raises(KeyError):
        sis.validate_input("CALC_PRESSURE", True)
    # Double
    assert sis.validate_input("ELEC_TEMP", "200")
    assert sis.validate_input("ELEC_TEMP", "200.0")
    assert sis.validate_input("ELEC_TEMP", 200.0)
    assert sis.validate_input("ELEC_TEMP", np.float64(200.0))
    # Integer array
    assert sis.validate_input("RELAX", "0 0 0")
    assert sis.validate_input("RELAX", "   0 0 0\n   0   0   0   ")
    assert sis.validate_input("RELAX", [0, 0, 0])
    # TODO: make sure this case doesn't throw warning
    assert sis.validate_input("RELAX", ["0 0 0", "0 0 0"])
    # Non-padded
    assert sis.validate_input("RELAX", [" 0    0  0", " 0 0   0 "])
    assert sis.validate_input("RELAX", np.array([0, 0, 0]))
    assert sis.validate_input("RELAX", np.array([True, False, True]))
    # TODO: exceptions for cases where int cannot be presented as bool!
    # Double array
    assert sis.validate_input("LATVEC", "1.0 0 0\n0 1.0 0\n0 0 1.0")
    assert sis.validate_input(
        "LATVEC", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    )
    assert sis.validate_input("LATVEC", np.eye(3))
    # TODO: extra user cases for character !


def test_api_validate_all_defaults():
    """All defaults given in the examples should be valid!"""
    from sparc.inputs import SparcInputs
    import numpy as np

    sis = SparcInputs()
    for param, pd in sis.parameters.items():
        default = pd.get("default", None)
        if default and (default != "auto"):
            assert sis.validate_input(param, default)


def test_api_convert_string():
    """Test if read string makes sense"""
    from sparc.inputs import SparcInputs
    import numpy as np

    sis = SparcInputs()
    # Integer
    assert sis.convert_string_to_value("CALC_PRES", "0") is False
    assert sis.convert_string_to_value("CALC_PRES", "0 ") == 0
    assert sis.convert_string_to_value("CALC_PRES", "1") is True
    # Double
    assert sis.convert_string_to_value("ELEC_TEMP", "200") == 200
    assert sis.convert_string_to_value("ELEC_TEMP", "200.0") == 200
    # Integer array
    relax_arr = sis.convert_string_to_value("RELAX", "0 0 0")
    print("relax,", relax_arr)
    assert len(relax_arr) == 3
    assert all([not p for p in relax_arr])
    relax_arr2 = sis.convert_string_to_value("RELAX", ["0 0 0", "0 0 0"])
    assert relax_arr2.ndim == 2
    # TODO: exceptions for cases where int cannot be presented as bool!
    # Double array
    latvec = sis.convert_string_to_value("LATVEC", "1.0 0 0\n0 1.0 0\n0 0 1.0")
    print("LATVEC", latvec)
    assert np.isclose(latvec, np.eye(3)).all()


def test_api_write_string():
    """Test if writing string from value makes sense"""
    from sparc.inputs import SparcInputs
    import numpy as np

    sis = SparcInputs()
    # Integer
    ref_s = "0"
    assert sis.convert_value_to_string("CALC_PRES", 0) == ref_s
    assert sis.convert_value_to_string("CALC_PRES", "0 ") == ref_s
    assert sis.convert_value_to_string("CALC_PRES", False) == ref_s
    # Double
    assert float(sis.convert_value_to_string("ELEC_TEMP", 200)) == 200
    assert float(sis.convert_value_to_string("ELEC_TEMP", 200.0)) == 200
    assert float(sis.convert_value_to_string("ELEC_TEMP", -200.0)) == -200
    # Integer array
    ref_s = "0 0 0"
    assert sis.convert_value_to_string("RELAX", [0, 0, 0]) == ref_s
    assert sis.convert_value_to_string("RELAX", [False, False, False]) == ref_s
    # TODO: exceptions for cases where int cannot be presented as bool!
    # Double array
    latvec_s = sis.convert_value_to_string("LATVEC", np.eye(3))
    assert latvec_s.count("\n") == 2
