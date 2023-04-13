import json
from pathlib import Path
from warnings import warn
import numpy as np

curdir = Path(__file__).parent
# TODO: must clean the api directory
default_api_dir = curdir / "sparc_json_api"
default_json_api = default_api_dir / "parameters.json"


class SparcInputs:
    def __init__(self, json_api=None):
        """Initialize the API from a json file"""
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
        # TODO: Make a parameters by categories

    def get_parameter_dict(self, parameter):
        parameter = parameter.upper()
        if parameter not in self.parameters.keys():
            raise KeyError(
                f"Parameter {parameter} is not known to SPARC {self.sparc_version}!"
            )
        return self.parameters[parameter]

    def help_info(self, parameter):
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
        """Give a string for a parameter, determine if the input follows the type

        input can be either a string or a 'direct' data type, like python float or numpy float

        TODO: there are many exceptions in array types, should enumerate
        """
        is_input_string = isinstance(input, str)

        pdict = self.get_parameter_dict(parameter)
        dtype = pdict["type"]
        if dtype == "string":
            return is_input_string
        elif dtype == "other":
            # Do nother for the "other" types but simply reply on the str() method
            if not is_input_string:
                warn(
                    f"Parameter {parameter} has 'other' data type and your input is not a string. I hope you know what you're doing!"
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
                return False
        elif "array" in dtype:
            if is_input_string:
                if ("." in input) and ("integer" in dtype):
                    warn(
                        (
                            f"Input {input} for parameter {parameter} it not strictly integer. "
                            "I can still perform the conversion but be aware of data loss"
                        )
                    )
                try:
                    arr = np.genfromtxt(input.splitlines(), delimiter=" ")
                except Exception:
                    arr = np.array(0.0)
            else:
                try:
                    arr = np.asarray(input)
                    if (arr.dtype != int) and ("integer" in dtype):
                        warn(
                            (
                                f"Input {input} for parameter {parameter} it not strictly intr. "
                                "I can still perform the conversion but be aware of data loss"
                            )
                        )
                except Exception:
                    arr = np.array(0.0)
            return len(arr.shape) > 0
        else:
            raise ValueError(f"Data type {dtype} is not supported!")

    def convert_string_to_value(self, parameter, string):
        """Convert a string input into valie parameter type"""
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
            value = float(string)
        elif dtype == "integer array":
            value = np.genfromtxt(string.splitlines(), delimiter=" ", dtype=int)
            if allow_bool_input:
                value = value.astype(bool)
        elif dtype == "double array":
            value = np.genfromtxt(string.splitlines(), delimiter=" ", dtype=float)
        else:
            # should not happen since validate_input has gatekeeping
            raise ValueError(f"Unsupported type {dtype}")

        return value
