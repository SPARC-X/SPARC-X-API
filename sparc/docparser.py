# -*- coding: utf-8 -*-
"""
A module to parse the latex documents provided by SPARC
and convert to its Python API

Created on Wed Mar  1 15:32:31 EST 2023

Tian Tian (alchem0x2a@gmail.com)
"""
import json
import re
from copy import copy
from datetime import datetime
from pathlib import Path
from warnings import warn

import numpy as np

# Some fields in master SPARC doc may cause auto type detection
# to fail, need hard-coded post-processing for now
postprocess_items = {
    "RELAX_FLAG": {"allow_bool_input": False},
    "NPT_SCALE_CONSTRAINTS": {"type": "string"},
    "NPT_SCALE_VECS": {"type": "integer array"},
    "TOL_POISSON": {"type": "double"},
}

sparc_repo_url = "https://github.com/SPARC-X/SPARC.git"


class SparcDocParser(object):
    """Parses LaTeX documentation of SPARC-X and converts it into a Python API.

    This class extracts parameter information from LaTeX source files,
    organizing it into a structured format that can be easily used in
    Python. It supports parsing of version details, parameter types,
    units, and other relevant information.

    Attributes:
        version (str): Parsed SPARC version, based on the documentation.
        parameter_categories (list): Categories of parameters extracted.
        parameters (dict): Extracted parameters with detailed information.
        other_parameters (dict): Additional parameters not categorized.

    Methods:
        find_main_file(main_file_pattern): Finds the main LaTeX file based on a pattern.
        get_include_files(): Retrieves a list of included LaTeX files.
        parse_version(parse): Parses and sets the SPARC version.
        parse_parameters(): Extracts parameters from LaTeX files.
        postprocess(): Applies hard-coded post-processing to some parameters.
        to_dict(): Converts parsed information into a dictionary.
        json_from_directory(directory, include_subdirs, **kwargs): Class method to create JSON from a directory.
        json_from_repo(url, version, include_subdirs, **kwargs): Class method to create JSON from a repository.

    """

    def __init__(
        self,
        directory=".",
        main_file="*Manual.tex",
        intro_file="Introduction.tex",
        params_from_intro=True,
        parse_version=True,
    ):
        """Create the doc parser pointing to the root of the doc file of SPARC

        The SPARC doc is organized as follows:
        SPARC/doc/.LaTeX/
            |---- Manual.tex
                  |---- Introduction.tex
                        |---- {Section}.tex

        For parameters additional to the standard SPARC options, such as the SQ / cyclix
        options, we merge the dict from the sub-dirs

        Args:
            doc_root: root directory to the LaTeX files, may look like `SPARC/doc/.LaTeX`
            main_file: main LaTeX file for the manual
            intro_file: LaTeX file for the introduction
            params_from_intro: only contain the parameters that can be parsed in `intro_file`
            parse_date: get the SPARC version by date
        """
        self.root = Path(directory)
        self.main_file = self.find_main_file(main_file)
        self.intro_file = self.root / intro_file
        if not self.intro_file.is_file():
            raise FileNotFoundError(f"Introduction file {intro_file} is missing!")
        self.include_files = self.get_include_files()
        self.params_from_intro = params_from_intro
        self.parse_version(parse_version)
        self.parse_parameters()
        self.postprocess()

    def find_main_file(self, main_file_pattern):
        """
        Finds the main LaTeX file that matches the given pattern, e.g. Manual.tex or Manual_cyclix.te

        Args:
            main_file_pattern (str): Pattern to match the main LaTeX file name.

        Returns:
            Path: Path to the main LaTeX file.

        Raises:
            FileNotFoundError: If no or multiple files match the pattern.
        """
        candidates = list(self.root.glob(main_file_pattern))
        if len(candidates) != 1:
            raise FileNotFoundError(
                f"Main file {main_file_pattern} is missing or more than 1 exists!"
            )
        return candidates[0]

    def get_include_files(self):
        """
        Retrieves a list of LaTeX files included in the main LaTeX document, e.g.  Manual.tex.

        Returns:
            list: A list of paths to the included LaTeX files.
        """
        pattern = r"\\begin\{document\}(.*?)\\end\{document\}"
        text = open(self.main_file, "r", encoding="utf8").read()
        # Only the first begin/end document will be matched
        match = re.findall(pattern, text, re.DOTALL)[0]
        pattern_include = r"\\include\{(.+?)\}"
        include = re.findall(pattern_include, match, re.DOTALL)
        include_files = []
        for name in include:
            tex_file = self.root / f"{name}.tex"
            if tex_file.is_file():
                include_files.append(tex_file)
            else:
                warn(
                    (
                        f"TeX file {tex_file} is missing! It may be a typo in the document, "
                        "ignore parameters from this file."
                    )
                )
        return include_files

    def parse_version(self, parse=True):
        """
        Parses and sets the SPARC version based on the C-source file, if possible.
        The date for the SPARC code is parsed from initialization.c in the "YYYY.MM.DD"
        format.

        Args:
            parse (bool): Whether to parse the version from the documentation.

        Sets:
            self.version (str): The parsed version in 'YYYY.MM.DD' format or None,
                                if either parse=False, or the C-source code is missing
        """
        if parse is False:
            self.version = None
            return
        init_c = self.root.parents[1] / "src" / "initialization.c"
        if not init_c.is_file():
            warn(
                'Cannot find the c source file "initialization.c", skip version parsing!'
            )
            self.version = None
            return
        text = open(init_c, "r", encoding="utf8").read()
        pattern_version = r"SPARC\s+\(\s*?version(.*?)\)"
        match = re.findall(pattern_version, text)
        if len(match) != 1:
            warn(
                'Parsing c source file "initialization.c" for version is unsuccessful!'
            )
            self.version = None
            return
        # We need to add more spacing matching in case the source code includes extra
        date_str = re.sub(r"\s+", " ", match[0].strip().replace(",", " "))
        # Older version of SPARC doc may contain abbreviated month format
        date_version = None
        for fmt in ("%b %d %Y", "%B %d %Y"):
            try:
                date_version = datetime.strptime(date_str, fmt).strftime("%Y.%m.%d")
                break
            except Exception:
                continue
        if date_version is None:
            raise ValueError(f"Cannot parse date time {date_str}")
        self.version = date_version
        return

    def __parse_parameter_from_frame(self, frame):
        """Parse the parameters from a single LaTeX frame

        Args:
            frame (str): a string containing the LaTeX frame (e.g. \\begin{frame} ... \\end{frame})

        Returns:
            dict: a key-value paired dict parsed from the frame. Some field names include:
                  name: TOL_POISSON
                  type: Double | Integer | String | Character | Double array
                  unit: specified in the doc
        """
        pattern_label = r"\\texttt\{(.*?)\}.*?\\label\{(.*?)\}"
        pattern_block = r"\\begin\{block\}\{(.*?)\}([\s\S]*?)\\end\{block\}"
        match_label = re.findall(pattern_label, frame, re.DOTALL | re.MULTILINE)
        if len(match_label) != 1:
            warn("Provided a non-structured frame for parsing, skip.")
            return {}
        symbol, label = (
            convert_tex_parameter(match_label[0][0].strip()),
            match_label[0][1].strip(),
        )
        # Every match contains the (name, content) pair of the blocks
        matches = re.findall(pattern_block, frame, re.DOTALL | re.MULTILINE)
        param_dict = {"symbol": symbol, "label": label}
        # TODO: add more type definition
        for key, content in matches:
            key = key.lower()
            content = content.strip()
            # Do not parse commented-out values

            if (key == "type") and (content.startswith("%")):
                warn(f"Parameter {symbol} is disabled in the doc, ignore!")
                return {}
            if key in ("example",):
                content = convert_tex_example(content)
            param_dict[key] = content
        # Sanitize 1: Convert types
        param_dict = sanitize_type(param_dict)
        # Sanitize 2: Convert default values
        param_dict = sanitize_default(param_dict)
        # Sanitize 3: Remove TeX components in description and remark
        param_dict = sanitize_description(param_dict)

        return param_dict

    def __parse_frames_from_text(self, text):
        """Extract all the frames that aren't commented in the text

        Arguments:
            text (str): Full LaTeX text
        Returns:
            list: Matched LaTeX Beamer frame fragments
        """
        pattern_frame = r"\\begin\{frame\}(.*?)\\end\{frame\}"
        matches = re.findall(pattern_frame, text, re.DOTALL | re.MULTILINE)
        return matches

    def __parse_intro_file(self):
        """Parse the introduction file

        Returns:
            parameter_dict (dict): dictionary using the parameter category as the main key
                            (following order in Introduction.tex)
            parameter_categories (list): list of categories
        """
        text_intro = open(self.intro_file, "r", encoding="utf8").read()
        pattern_params = (
            r"^\\begin\{frame\}.*?\{Input file options\}.*?$(.*?)\\end\{frame\}"
        )
        pattern_block = r"\\begin\{block\}\{(.*?)\}([\s\S]*?)\\end\{block\}"
        pattern_line = r"\\hyperlink\{(.*?)\}{\\texttt\{(.*?)\}\}"
        text_params = re.findall(pattern_params, text_intro, re.DOTALL | re.MULTILINE)[
            0
        ]
        parameter_categories = []
        parameter_dict = {}
        for match in re.findall(pattern_block, text_params):
            cat = match[0].lower()
            # print(cat)
            if cat in parameter_categories:
                raise ValueError(
                    f"Key {cat} already exists! You might have a wrong LaTeX doc file!"
                )
            parameter_categories.append(cat)
            parameter_dict[cat] = []
            param_lines = match[1].split("\n")
            for line in param_lines:
                matches = re.findall(pattern_line, line)
                if len(matches) == 0:
                    continue
                # Each match should contain 2 items, the "Link" that matches a reference in included-tex files
                # symbol is the actual symbol name (in text-format)
                # In most cases the link and symbol should be the same
                for match in matches:
                    label, symbol = match[0].strip(), convert_tex_parameter(
                        match[1].strip()
                    )
                    parameter_dict[cat].append({"label": label, "symbol": symbol})
        return parameter_categories, parameter_dict

    def __parse_all_included_files(self):
        """Pop up all known parameters from included files
        Returns:
            dict: All known parameters from included files
        """
        all_params = {}
        for f in self.include_files:
            # Do not parse intro file since it's waste of time
            if f.resolve() == self.intro_file.resolve():
                continue
            text = open(f, "r", encoding="utf8").read()
            frames = self.__parse_frames_from_text(text)
            for frame in frames:
                dic = self.__parse_parameter_from_frame(frame)
                if len(dic) > 0:
                    label = dic["label"]
                    all_params[label] = dic
        return all_params

    def parse_parameters(self):
        """The actual thing for parsing parameters

        Sets:
            parameters (dict): All parsed parameters
            parameter_categoris (list): List of categories
            other_parameters (dict): Any parameters that are not included in the categories
        """
        parameter_categories, parameter_dict = self.__parse_intro_file()
        all_params = self.__parse_all_included_files()
        self.parameter_categories = parameter_categories
        # parameters contain only the "valid" ones that are shown in the intro
        # all others are clustered in "other_parameters"
        self.parameters = {}
        for cat, params in parameter_dict.items():
            for p in params:
                label = p["label"]
                symbol = p["symbol"]
                param_details = all_params.pop(label, {})
                if param_details != {}:
                    param_details["category"] = cat
                    self.parameters[symbol] = param_details

        self.other_parameters = {}
        for param_details in all_params.values():
            symbol = param_details["symbol"]
            self.other_parameters[symbol] = param_details
        return

    def postprocess(self):
        """Use the hardcoded dict prostprocess_items to fix some issues"""
        for param, fix in postprocess_items.items():
            if param in self.parameters:
                self.parameters[param].update(**fix)
        return

    def to_dict(self):
        """Output a json dict from current document parser

        Returns:
            dict: All API schemes in dict
        """
        doc = {}
        doc["sparc_version"] = self.version
        doc["categories"] = self.parameter_categories
        doc["parameters"] = {k: v for k, v in sorted(self.parameters.items())}
        doc["other_parameters"] = {
            k: v for k, v in sorted(self.other_parameters.items())
        }
        doc["data_types"] = sorted(set([p["type"] for p in self.parameters.values()]))
        return doc

    @classmethod
    def json_from_directory(cls, directory=".", include_subdirs=True, **kwargs):
        """
        Recursively add parameters from all Manual files
        Arguments:
            directory (str or PosixPath): The directory to the LaTeX files, e.g. <sparc-root>/doc/.LaTeX
            include_subdirs (bool): If true, also parse the manual files in submodules, e.g. cyclix, highT
        Returns:
            str: Formatted json-string of the API
        """
        directory = Path(directory)
        root_dict = cls(directory=directory, **kwargs).to_dict()
        if include_subdirs:
            for sub_manual_tex in directory.glob("*/*Manual.tex"):
                subdir = sub_manual_tex.parent
                try:
                    sub_dict = cls(directory=subdir, parse_version=False).to_dict()
                except FileNotFoundError:
                    print(subdir, " Latex files not found. Check naming conventions for Manual.tex. Expects format *Manual.tex")
                    continue
                for param, param_desc in sub_dict["parameters"].items():
                    if param not in root_dict["parameters"]:
                        root_dict["parameters"][param] = param_desc
        json_string = json.dumps(root_dict, indent=True)
        return json_string

    @classmethod
    def json_from_repo(
        cls, url=sparc_repo_url, version="master", include_subdirs=True, **kwargs
    ):
        """
        Download the source code from git and use json_from_directory to parse
        Arguments:
            url (str): URL for the repository of SPARC, default is "https://github.com/SPARC-X/SPARC.git"
            version (str): Git version or commit hash of the SPARC repo
            include_subdirs (bool): If true, also parse the manual files in submodules, e.g. cyclix, highT
        Returns:
            str: Formatted json-string of the API
        """
        import tempfile
        from subprocess import run

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            download_dir = tmpdir / "SPARC"
            download_cmds = ["git", "clone", "--depth", "1", str(url), "SPARC"]
            run(download_cmds, cwd=tmpdir)
            if version not in ["master", "HEAD"]:
                fetch_cmds = ["git", "fetch", "--depth", "1", str(version)]
                run(fetch_cmds, cwd=download_dir)
                checkout_cmds = ["git", "checkout", str(version)]
                run(checkout_cmds, cwd=download_dir)
            json_string = cls.json_from_directory(
                directory=download_dir / "doc" / ".LaTeX",
                include_subdirs=include_subdirs,
                **kwargs,
            )
        return json_string


def convert_tex_parameter(text):
    """
    Conver a TeX string to non-escaped name (for parameter only)
    Arguments:
        text (str): Parameter name in LaTeX format
    Returns:
        str: Text with sanitized parameter
    """
    return text.strip().replace("\_", "_")


def convert_tex_example(text):
    """Convert TeX codes of examples as much as possible
    The examples follow the format
    SYMBOL: values (may contain new lines)
    Arguments:
        text (str): Single or multiline LaTeX contents
    Returns:
        str: Sanitized literal text
    """
    mapper = {"\\texttt{": "", "\_": "_", "}": "", "\\": "\n"}
    new_text = copy(text)
    for m, r in mapper.items():
        new_text = new_text.replace(m, r)

    symbol, values = new_text.split(":")
    symbol = symbol.strip()
    values = re.sub("\n+", "\n", values.strip())
    # Remove all comment lines
    values = "\n".join(
        [l for l in values.splitlines() if not l.lstrip().startswith("%")]
    )
    new_text = f"{symbol}: {values}"
    return new_text


def convert_tex_default(text, desired_type=None):
    """Convert default values as much as possible.
    The desire type will convert the default values
    to the closest format

    Currently supported conversions
    1. Remove all surrounding text modifiers (texttt)
    2. Remove all symbol wrappers $
    3. Convert value to single or array

    Arguments:
        text (str): Raw text string for value
        desired_type (str or None): Data type to be converted to. If None, preserve the string format

    Returns:
        converted: Value converted from raw text
    """
    mapper = {
        "\\texttt{": "",
        "}": "",
        "{": "",
        "\\_": "_",
        "\_": "_",
        "\\\\": "\n",
        "$": "",
    }
    text = text.strip()
    text = re.sub(r"\\hyperlink\{.*?\}", "", text)
    text = re.sub(r"\\times", "x", text)
    for m, r in mapper.items():
        text = text.replace(m, r)
    text = re.sub(r"\n+", "\n", text)
    # Remove all comment lines
    text = "\n".join([l for l in text.splitlines() if not l.lstrip().startswith("%")])

    # print(text)
    converted = None
    if "none" in text.lower():
        converted = None
    elif "no default" in text.lower():
        converted = None
    elif "automat" in text.lower():
        converted = "auto"
    else:
        # try type conversion
        if desired_type is None:
            converted = text
        elif desired_type == "string":
            converted = text
        else:
            converted = text2value(text, desired_type)
    return converted


def convert_comment(text):
    """Used to remove TeX-specific commands in description and remarks
    as much as possible

    Arguments:
        text (str): Raw LaTeX code for the comment section in manual

    Returns:
        str: Sanitized plain text
    """
    mapper = {
        "\\texttt{": "",
        "}": "",
        "{": "",
        "\\_": "_",
        "\_": "_",
        "\\\\": "\n",
        "$": "",
    }
    text = text.strip()
    text = re.sub(r"\\hyperlink\{.*?\}", "", text)
    text = re.sub(r"\\href\{.*?\}", "", text)
    text = re.sub(r"\\times", "x", text)
    for m, r in mapper.items():
        text = text.replace(m, r)
    text = re.sub(r"\n+", "\n", text)
    # Remove all comment lines
    text = "\n".join([l for l in text.splitlines() if not l.lstrip().startswith("%")])
    return text


def text2value(text, desired_type):
    """Convert raw text to a desired type

    Arguments:
        text (str): Text contents for the value
        desired_type (str): Target data type from 'string', 'integer',
                            'integer array', 'double', 'double array',
                            'bool', 'bool array'
    Returns:
        converted: Value converted to the desired type
    """
    if desired_type is None:
        return text
    desired_type = desired_type.lower()
    if desired_type == "string":
        return text.strip()

    try:
        arr = np.genfromtxt(text.splitlines(), delimiter=" ", dtype=float)
        if np.isnan(arr).any():
            warn(
                f"Some fields in {text} cannot converted to a numerical array, will skip conversion."
            )
            arr = None
    except Exception as e:
        warn(
            f"Cannot transform {text} to array, skip converting. Error message is:\n {e}"
        )
        arr = None

    if arr is None:
        return None

    # Upshape ndarray to at least 1D
    if arr.shape == ():
        arr = np.reshape(arr, [1])

    converted = None
    from contextlib import suppress

    # Ignore all failures and make conversion None
    with suppress(Exception):
        if desired_type == "integer":
            converted = int(arr[0])
        elif desired_type == "bool":
            converted = bool(arr[0])
        elif desired_type == "double":
            converted = float(arr[0])
        elif desired_type == "integer array":
            converted = np.ndarray.tolist(arr.astype(int))
        elif desired_type == "bool array":
            converted = np.ndarray.tolist(arr.astype(bool))
        elif desired_type == "double array":
            converted = np.ndarray.tolist(arr.astype(float))
    return converted


def is_array(text):
    """Simply try to convert a string into a numpy array and compare if length is larger than 1
    it is only used to compare a float / int value
    """
    val = np.fromstring(text, sep=" ")
    if len(val) == 1:
        return False
    else:
        return True


def contain_only_bool(text):
    """Check if a string only contains 0 1 or spaces"""
    if any([c in text for c in (".", "+", "-", "e", "E")]):
        return False
    digits = re.findall(r"[-+e\d]+", text, re.DOTALL)
    for d in digits:
        val = int(d)
        if val not in (0, 1):
            return False
    return True


def sanitize_description(param_dict):
    """Sanitize the description and remark field

    Arguments:
        param_dict (dict): Raw dict for one parameter entry

    Returns:
        dict: Sanitized parameter dict with comment, remark and description
              converted to human-readable formats
    """
    sanitized_dict = param_dict.copy()

    original_desc = sanitized_dict["description"]
    sanitized_dict["description_raw"] = original_desc

    original_remark = sanitized_dict.get("remark", "")
    sanitized_dict["remark_raw"] = original_remark

    sanitized_dict["description"] = convert_comment(original_desc)
    sanitized_dict["remark"] = convert_comment(original_remark)
    return sanitized_dict


def sanitize_default(param_dict):
    """Sanitize the default field
    1. Create an extra field `default_remark` that copies original default
    2. Use `convert_tex_default` to convert values as much as possible

    This function should be called after sanitize_type
    """
    sanitized_dict = param_dict.copy()
    original_default = sanitized_dict["default"]
    sanitized_dict["default_remark"] = original_default
    converted_default = convert_tex_default(original_default, param_dict["type"])
    sanitized_dict["default"] = converted_default
    return sanitized_dict


def sanitize_type(param_dict):
    """Sanitize the param dict so that the type are more consistent

    For example, if type is Double / Integer,
    but parameter is a vector,
    make a double vector or integer vector
    """
    sanitized_dict = param_dict.copy()
    symbol = param_dict["symbol"]
    origin_type = param_dict.get("type", None)
    if origin_type is None:
        print("Dict does not have type!")
        return sanitized_dict
    origin_type = origin_type.lower()

    sanitized_type = None
    sanitized_dict["allow_bool_input"] = False
    # First pass, remove all singular types
    if origin_type == "0 or 1":
        origin_type = "integer"
    elif "permutation" in origin_type:
        sanitized_type = "integer"
    elif origin_type in ("string", "character"):
        sanitized_type = "string"
    elif "array" in origin_type:
        sanitized_type = origin_type

    # Pass 2, test if int values are arrays
    if (origin_type in ["int", "integer", "double"]) and (sanitized_type is None):
        if "int" in origin_type:
            origin_type = "integer"
        # Test if the value from example is a single value or array
        try:
            example_value = param_dict["example"].split(":")[1]
            default = param_dict["default"]
            _array_test = is_array(example_value)
            _bool_test = contain_only_bool(example_value) and contain_only_bool(default)
        except Exception as e:
            warn(
                f"Array conversion failed for {example_value}, ignore."
                f"The error is {e}"
            )
            _array_test = False  # Retain

        if _array_test is True:
            sanitized_type = f"{origin_type} array"
        else:
            sanitized_type = origin_type

        # Pass 3: int to boolean test. This should be done very tight
        if _bool_test and ("integer" in sanitized_type):
            sanitized_dict["allow_bool_input"] = True

    if sanitized_type is None:
        # Currently there is only one NPT_NH_QMASS has this type
        # TODO: think of a way to format a mixed array?
        warn(f"Type of {symbol} if not standard digit or array, mark as others.")
        sanitized_type = "other"
        # TODO: how about provide a true / false type?
    sanitized_dict["type"] = sanitized_type
    return sanitized_dict


if __name__ == "__main__":
    # Run the module as independent script to extract a json-formatted parameter list
    from argparse import ArgumentParser

    argp = ArgumentParser(description="Parse the LaTeX doc to json")
    argp.add_argument(
        "-o",
        "--output",
        default="parameters.json",
        help="Output file name (json-formatted)",
    )
    argp.add_argument(
        "--include-subdirs",
        action="store_true",
        help="Parse manual parameters from subdirs",
    )
    argp.add_argument("--git", action="store_true")
    argp.add_argument(
        "--version",
        default="master",
        help="Version of the doc. Only works when using git repo",
    )
    argp.add_argument(
        "root",
        nargs="?",
        help=(
            "Root of the SPARC doc LaTeX files, or remote git repo link. If not provided and --git is enables, use the default github repo"
        ),
    )

    args = argp.parse_args()
    output = Path(args.output).with_suffix(".json")
    if args.git:
        if args.root is None:
            root = sparc_repo_url
        else:
            root = args.root
        json_string = SparcDocParser.json_from_repo(
            url=root, version=args.version, include_subdirs=args.include_subdirs
        )
    else:
        json_string = SparcDocParser.json_from_directory(
            directory=Path(args.root), include_subdirs=args.include_subdirs
        )
    with open(output, "w", encoding="utf8") as fd:
        fd.write(json_string)
    print(f"SPARC parameter specifications written to {output}!")
    print("If you need to fintune the definitions, please edit them manually.")
