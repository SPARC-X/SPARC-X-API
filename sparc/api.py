import json
from io import StringIO
from pathlib import Path
from warnings import warn

import numpy as np

curdir = Path(__file__).parent
default_api_dir = curdir / "sparc_json_api"
default_json_api = default_api_dir / "parameters.json"


class SparcAPI:
    """
    An interface to the parameter settings in SPARC-X calculator. User can use the
    SparcAPI instance to validate and translate parameters that matches a certain
    version of the SPARC-X code.

    Attributes:
        sparc_version (str): Version of SPARC.
        categories (dict): Categories of parameters.
        parameters (dict): Detailed parameters information.
        other_parameters (dict): Additional parameters.
        data_types (dict): Supported data types.

    Methods:
        get_parameter_dict(parameter): Retrieves dictionary for a specific parameter.
        help_info(parameter): Provides detailed information about a parameter.
        validate_input(parameter, input): Validates user input against the expected parameter type.
        convert_string_to_value(parameter, string): Converts string input to the appropriate data type.
        convert_value_to_string(parameter, value): Converts a value to a string representation.
    """

    def __init__(self, json_api=None):
        """ """
        if json_api is None:
            json_api = Path(default_json_api)
        else:
            json_api = Path(json_api)

        json_data = json.load(open(json_api, "r"))
        self.sparc_version = json_data["sparc_version"]
        self.categories = json_data["categories"]
        self.parameters = json_data["parameters"]
        self.other_parameters = json_data["other_parameters"]
        self.data_types = json_data["data_types"]
        # TT: 2024-10-31 add the sources to trace the origin
        # locate_api can modify self.source if it is deferred from LaTeX
        # at runtime
        self.source = {"path": json_api.as_posix(), "type": "json"}

    def get_parameter_dict(self, parameter):
        """
        Retrieves the dictionary for a specified parameter.

        Args:
            parameter (str): The name of the parameter.

        Returns:
            dict: Dictionary containing details of the parameter.

        Raises:
            KeyError: If the parameter is not known to the SPARC version.
        """
        parameter = parameter.upper()
        if parameter not in self.parameters.keys():
            raise KeyError(
                f"Parameter {parameter} is not known to " f"SPARC {self.sparc_version}!"
            )
        return self.parameters[parameter]

    def help_info(self, parameter):
        """Provides a detailed information string for a given parameter.

        Args:
            parameter (str): The name of the parameter to get information for.

        Returns:
            str: A formatted string with detailed information about the parameter.
        """
        pdict = self.get_parameter_dict(parameter)
        message = "\n".join(
            [
                f"{key}: {pdict[key]}"
                for key in (
                    "symbol",
                    "category",
                    "type",
                    "unit",
                    "default",
                    "example",
                    "description",
                    "remark",
                    "allow_bool_input",
                )
            ]
        )
        return message

    def validate_input(self, parameter, input):
        """
        Validates if the given input is appropriate for the specified parameter's type.

        Args:
            parameter (str): The name of the parameter.
            input: The input to validate, can be of various types (string, int, float, numpy types).

        Returns:
            bool: True if input is valid, False otherwise.

        Raises:
            ValueError: If the data type of the parameter is not supported.
        """
        is_input_string = isinstance(input, str)
        pdict = self.get_parameter_dict(parameter)
        dtype = pdict["type"]
        if dtype == "string":
            return is_input_string
        elif dtype == "other":
            # Do nother for the "other" types but simply
            # reply on the str() method
            if not is_input_string:
                warn(
                    f"Parameter {parameter} has 'other' data type "
                    "and your input is not a string. "
                    "I hope you know what you're doing!"
                )
            return True
        elif dtype == "integer":
            try:
                int(input)
                return True
            except (TypeError, ValueError):
                return False
        elif dtype == "double":
            try:
                float(input)
                return True
            except (TypeError, ValueError):
                try:
                    float(input.split()[0])
                    return True
                except Exception:
                    return False
        elif "array" in dtype:
            if is_input_string:
                if ("." in input) and ("integer" in dtype):
                    warn(
                        (
                            f"Input {input} for parameter "
                            f"{parameter} it not strictly integer. "
                            "I may still perform the conversion "
                            "but be aware of data loss"
                        )
                    )
                try:
                    arr = np.genfromtxt(input.splitlines(), dtype=float, ndmin=1)
                    # In valid input with nan
                    if np.isnan(arr).any():
                        arr = np.array(0.0)
                except Exception:
                    arr = np.array(0.0)
            else:
                try:
                    arr = np.atleast_1d(np.asarray(input))
                    if (arr.dtype not in (int, bool)) and ("integer" in dtype):
                        warn(
                            (
                                f"Input {input} for parameter {parameter} is"
                                " not strictly integer. "
                                "I can still perform the conversion but "
                                "be aware of data loss"
                            )
                        )
                except Exception:
                    arr = np.array(0.0)
            return len(arr.shape) > 0
        else:
            raise ValueError(f"Data type {dtype} is not supported!")

    def convert_string_to_value(self, parameter, string):
        """
        Converts a string input to the appropriate value type of the parameter.

        Args:
            parameter (str): The name of the parameter.
            string (str): The string input to convert.

        Returns:
            The converted value, type depends on parameter's expected type.

        Raises:
            TypeError: If the input is not a string.
            ValueError: If the string is not a valid input for the parameter.
        """

        # Special case, the string may be a multiline string-array!
        if isinstance(string, list):
            # Make sure there is a line break at the end, for cases like ["2."]
            string.append("")
            string = [s.strip() for s in string]
            string = "\n".join(string)

        is_input_string = isinstance(string, str)
        if not is_input_string:
            raise TypeError("Please give a string input!")

        if not self.validate_input(parameter, string):
            raise ValueError(f"{string} is not a valid input for {parameter}")

        pdict = self.get_parameter_dict(parameter)
        dtype = pdict["type"]
        allow_bool_input = pdict.get("allow_bool_input", False)

        if dtype == "string":
            value = string.strip()
        elif dtype == "integer":
            value = int(string)
            if allow_bool_input:
                value = bool(value)
        elif dtype == "double":
            # Some inputs, like TARGET_PRESSURE, may be accepted with a unit
            # like 0.0 GPa. Only accept the first part
            try:
                value = float(string)
            except ValueError as e:
                try:
                    value = float(string.split()[0])
                except Exception:
                    raise e
        elif dtype == "integer array":
            value = np.genfromtxt(string.splitlines(), dtype=int, ndmin=1)
            if allow_bool_input:
                value = value.astype(bool)
        elif dtype == "double array":
            value = np.genfromtxt(string.splitlines(), dtype=float, ndmin=1)
        elif dtype == "other":
            value = string
            # should not happen since validate_input has gatekeeping
        else:
            raise ValueError(f"Unsupported type {dtype}")

        return value

    def convert_value_to_string(self, parameter, value):
        """
        Converts a value to its string representation based on the parameter type.

        Args:
            parameter (str): The name of the parameter.
            value: The value to convert.

        Returns:
            str: The string representation of the value.

        Raises:
            ValueError: If the value is not valid for the parameter.
        """

        is_input_string = isinstance(value, str)
        if not self.validate_input(parameter, value):
            raise ValueError(f"{value} is not a valid input for {parameter}")

        # Do not conver, just return the non-padded string
        if is_input_string:
            return value.strip()

        pdict = self.get_parameter_dict(parameter)
        dtype = pdict["type"]
        # allow_bool_input = pdict.get("allow_bool_input", False)

        if dtype == "string":
            string = str(value).strip()
        elif dtype == "integer":
            # Be aware of bool values!
            string = str(int(value))
        elif dtype == "double":
            string = "{:.14f}".format(float(value))
        elif dtype in ("integer array", "double array"):
            string = _array_to_string(value, dtype)
        elif dtype == "other":
            if not is_input_string:
                raise ValueError("Only support string value when datatype is other")
            string = value
        else:
            # should not happen since validate_input has gatekeeping
            raise ValueError(f"Unsupported type {dtype}")

        return string


def _array_to_string(arr, format):
    """
    Converts an array to a string representation based on the specified format.

    Args:
        arr (array): The array to convert.
        format (str): The format type ('integer array', 'double array', etc.).

    Returns:
        str: String representation of the array.
    """
    arr = np.array(arr)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    buf = StringIO()
    if format in ("integer array", "integer"):
        fmt = "%d"
    elif format in ("double array", "double"):
        fmt = "%.14f"
    np.savetxt(buf, arr, delimiter=" ", fmt=fmt, header="", footer="", newline="\n")
    # Return the string output of the buffer with
    # whitespaces removed
    return buf.getvalue().strip()
