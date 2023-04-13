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
    # Double
    assert sis.validate_input("ELEC_TEMP", "200")
    assert sis.validate_input("ELEC_TEMP", "200.0")
    assert sis.validate_input("ELEC_TEMP", 200.0)
    assert sis.validate_input("ELEC_TEMP", np.float64(200.0))
    # Integer array
    assert sis.validate_input("RELAX", "0 0 0")
    assert sis.validate_input("RELAX", [0, 0, 0])
    assert sis.validate_input("RELAX", np.array([0, 0, 0]))
    assert sis.validate_input("RELAX", np.array([True, False, True]))
    # TODO: exceptions for cases where int cannot be presented as bool!
    # Double array
    assert sis.validate_input("LATVEC", "1.0 0 0\n0 1.0 0\n0 0 1.0")
    assert sis.validate_input("LATVEC", [[1.0, 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    assert sis.validate_input("LATVEC", np.eye(3))
    # TODO: extra user cases for character !

def test_api_validate_all_defaults():
    """All defaults given in the examples should be valid!
    """
    from sparc.inputs import SparcInputs
    import numpy as np
    sis = SparcInputs()
    for param, pd in sis.parameters.items():
        default = pd.get("default", None)
        if default and (default != "auto"):
            assert sis.validate_input(param, default)
