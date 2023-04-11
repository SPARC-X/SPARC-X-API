# -*- coding: utf-8 -*-
"""
A module to parse the latex documents provided by SPARC
and convert to its Python API

Created on Wed Mar  1 15:32:31 EST 2023

Tian Tian (alchem0x2a@gmail.com)
"""
import re
import os
import json
from pathlib import Path
from warnings import warn
import numpy as np
from copy import copy


class SPARCDocParser(object):
    """Use regex to parse LaTeX doc to python API
    """

    def __init__(self, doc_root=".",
                 main_file="Manual.tex",
                 intro_file="Introduction.tex",
                 params_from_intro=True):
        """Create the doc parser pointing to the root of the doc file of SPARC

        The SPARC doc is organized as follows:
        SPARC/doc/.LaTeX/
            |---- Manual.tex
                  |---- Introduction.tex  
                        |---- {Section}.tex
        TODO: include the parameters for SQ / HT calculations

        Arguments:
        `doc_root`: root directory to the LaTeX files, may look like `SPARC/doc/.LaTeX`
        `main_file`: main LaTeX file for the manual
        `intro_file`: LaTeX file for the introduction
        `params_from_intro`: only contain the parameters that can be parsed in `intro_file`
        """
        self.root = Path(doc_root)
        self.main_file = self.root / main_file
        if not self.main_file.is_file():
            raise FileNotFoundError(f"Main file {main_file} is missing!")
        self.intro_file = self.root / intro_file
        if not self.intro_file.is_file():
            raise FileNotFoundError(f"Introduction file {intro_file} is missing!")
        self.include_files = self.get_include_files()
        self.params_from_intro = params_from_intro
        self.param_dict = self.parse_parameters()

    def get_include_files(self):
        """Get a list of included LaTeX files from Manual.tex
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
                warn((f"TeX file {tex_file} is missing! It may be a typo in the document, "
                      "ignore parameters from this file."))
        return include_files

    def __parse_parameter_from_frame(self, frame):
        """Parse the parameters from a single LaTeX frame

        Arguments:
        `frame`: a string containing the LaTeX frame (e.g. \begin{frame} ... \end{frame})
                 
        fields are:
        name: TOL_POISSON
        type: Double | Integer | String | Character | Double array
        unit: specified in the doc
        """
        pattern_label = r"\\texttt\{(.*?)\}.*?\\label\{(.*?)\}"
        pattern_block = r"\\begin\{block\}\{(.*?)\}([\s\S]*?)\\end\{block\}"
        match_label = re.findall(pattern_label, frame, re.DOTALL | re.MULTILINE)
        if len(match_label) != 1:
            # warn("Provided a non-structured frame for parsing, skip.")
            return {}
        # print(match_label)
        symbol, label = convert_tex_parameter(match_label[0][0].strip()), match_label[0][1].strip()
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
        param_dict = sanitize_type(param_dict)
        return param_dict

    def __parse_frames_from_text(self, text):
        """Extract all the frames that aren't commented in the text

        Arguments:
        `text`: LaTeX text
        """
        pattern_frame = r"\\begin\{frame\}(.*?)\\end\{frame\}"
        matches = re.findall(pattern_frame, text, re.DOTALL | re.MULTILINE)
        return matches
    
    def __parse_intro_file(self):
        """Parse the introduction file

        Returns:
        `parameter_dict`: dictionary using the parameter category as the main key 
                          (following order in Introduction.tex)
        `parameter_categories`: list of categories
        """
        text_intro = open(self.intro_file, "r", encoding="utf8").read()
        # import pdb; pdb.set_trace()
        pattern_params = r"^\\begin\{frame\}.*?\{Input file options\}.*?$(.*?)\\end\{frame\}"
        pattern_block =  r"\\begin\{block\}\{(.*?)\}([\s\S]*?)\\end\{block\}"
        pattern_line = r"\\hyperlink\{(.*?)\}{\\texttt\{(.*?)\}\}"
        text_params = re.findall(pattern_params, text_intro, re.DOTALL | re.MULTILINE)[0]
        parameter_categories = []
        parameter_dict = {}
        for match in re.findall(pattern_block, text_params):
            cat = match[0].lower()
            print(cat)
            if cat in parameter_categories:
                raise ValueError(f"Key {cat} already exists! You might have a wrong LaTeX doc file!")
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
                    label, symbol = match[0].strip(), convert_tex_parameter(match[1].strip())
                    print(label, symbol)
                    parameter_dict[cat].append((label, symbol))
        return parameter_categories, parameter_dict

    def __parse_all_included_files(self):
        """Pop up all known parameters from included files,
        """
        all_params = {}
        for f in self.include_files:
            # Do not parse intro file since it's waste of time
            if f.resolve() == self.intro_file.resolve():
                continue
            print("Parsing", f)
            text = open(f, "r", encoding="utf8").read()
            frames = self.__parse_frames_from_text(text)
            # print(frames)
            for frame in frames:
                dic = self.__parse_parameter_from_frame(frame)
                # print(frame)
                if len(dic) > 0:
                    label = dic["label"]
                    all_params[label] = dic
        return all_params
                # print(dic)
        

    
    def parse_parameters(self):
        """The actual thing for parsing parameters
        """
        parameter_categories, parameter_dict = self.__parse_intro_file()
        all_params = self.__parse_all_included_files()
        

def convert_tex_parameter(text):
    """Conver a TeX string to non-escaped name (for parameter only)
    """
    return text.strip().replace("\_", "_")

def convert_tex_example(text):
    """Convert TeX codes of examples as much as possible
    The examples follow the format
    SYMBOL: values (may contain new lines)
    """
    mapper = {"\\texttt{": "", "\_": "_", "}": "", "\\": "\n"}
    new_text = copy(text)
    for m, r in mapper.items():
        new_text = new_text.replace(m, r)

    symbol, values = new_text.split(":")
    symbol = symbol.strip()
    values = re.sub('\n+', '\n', values.strip())
    new_text = f"{symbol}: {values}"
    return new_text

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
    """Check if a string only contains 0 1 or spaces
    """
    if any([c in text for c in (".", "+", "-", "e", "E")]):
        return False
    digits = re.findall(r"[-+e\d]+", text, re.DOTALL)
    for d in digits:
        val = int(d)
        if val not in (0, 1):
            return False
    return True
    

def sanitize_type(param_dict):
    """Sanitize the param dict so that the type are more consistent

    For example, if type is Double / Integer, but parameter is a vector, make a double vector or integer vector
    """
    sanitized_dict = param_dict.copy()
    symbol = param_dict["symbol"]
    origin_type = param_dict.get("type", None)
    if origin_type is None:
        print(f"Dict does not have type!")
        return sanitized_dict
    origin_type = origin_type.lower()

    sanitized_type = None
    # import pdb; pdb.set_trace()
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
            print("Example", param_dict["example"], example_value)
            # print()
            _array_test = is_array(example_value)
            _bool_test = contain_only_bool(example_value)
            print(_array_test)
        except Exception:
            raise
            _array_test = False # Retain
            
        if _array_test is True:
            sanitized_type = f"{origin_type} array"
        else:
            sanitized_type = origin_type

        # Pass 3: int to boolean test. Do not convert double!
        if _bool_test and ("integer" in sanitized_type):
            sanitized_type = sanitized_type.replace("integer", "bool")

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
    argp.add_argument("-o", "--output", default="parameters.json", help="Output file name (json-formatted)")
    argp.add_argument("root", help="Root directory of the latex files")  # root directory of the LaTeX files
    args = argp.parse_args()
    parser = SPARCDocParser(Path(args.root))
    print("The following files are included in the introduction:")
    for file in parser.include_files:
        print(file, file.exists())
