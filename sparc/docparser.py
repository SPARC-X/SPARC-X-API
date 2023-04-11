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
        pattern_block = r"\\begin\{block\}\{(.*?)\}([\s\S]*?)\\end\{block\}"
        # Every match contains the (name, content) pair of the blocks
        matches = re.findall(pattern_block, frame, re.DOTALL)
        param_dict = {}
        for key, content in matches:
            param_dict[key.lower()] = content
        return param_dict

    def __parse_frames_from_text(self, text):
        """Extract all the frames that aren't commented in the text

        Arguments:
        `text`: LaTeX text
        """
        pattern_frame = r"\\begin\{frame\}(?:\[(.*?))\])?\{(.*?)\}\\end\{frame\}"
        matches = re.findall(pattern_frame, text, re.DOTALL)
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
                    link, symbol = match[0].strip(), convert_tex_parameter(match[1].strip())
                    print(link, symbol)
                    parameter_dict[cat].append((link, symbol))
        return parameter_categories, parameter_dict

    
    def parse_parameters(self):
        """The actual thing for parsing parameters
        """
        parameter_categories, parameter_dict = self.__parse_intro_file()
        

def convert_tex_parameter(text):
    """Conver a TeX string to non-escaped name (for parameter only)
    """
    return text.strip().replace("\_", "_")

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
