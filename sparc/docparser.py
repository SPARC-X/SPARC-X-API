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

class SPARCDocParser(object):
    """Use regex to parse LaTeX doc to python API
    """
    def __init__(self, doc_root=".",
                main_file="Manual.tex",
                intro_file="Introduction.tex",
                params_from_intro=True):
        """Create the doc parser pointing to the root of the doc file

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
        self.include_files = self.get_include_files()
        

    def get_include_files(self):
        """Get a list of included LaTeX files from Manual.tex
        """
        pattern = r"\\begin\{document\}(.*?)\\end\{document\}"
        text = open(self.main_file, "r").read()
        # Only the first begin/end document will be matched
        match = re.findall(pattern, text, re.DOTALL)[0]
        pattern_include = r"\\include\{(.+?)\}"
        include = re.findall(pattern_include, match, re.DOTALL)
        include_files = [self.root / f"{name}.tex" for name in include]
        return include_files


    def __parse_parameter_from_frame(self, frame):
        """Parse the parameters from a single LaTeX frame
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
        """
        pattern_frame = r"\\begin\{frame\}(?:\[(.*?))\])?\{(.*?)\}\\end\{frame\}"
        matches = re.findall(pattern_frame, text, re.DOTALL)
        return matches
        
        

        
        
        
