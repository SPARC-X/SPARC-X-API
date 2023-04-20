"""Test for SPARC doc parser
"""
import pytest
from pathlib import Path
import os

curdir = Path(__file__).parent
test_doc_dir = curdir / "sparc-latex-doc-202302"


def test_docparser_init_wrong(fs):
    """Mimic situations where docparser inits at wrong file structure
    """
    from sparc.docparser import SPARCDocParser

    # Case 1: not a doc structure
    with pytest.raises(FileNotFoundError):
        sp = SPARCDocParser("/tmp")
    
    # Case 2: missing manual
    fs.create_dir(".LaTeX")
    fs.create_file(".LaTeX/Manual.tex")
    with pytest.raises(FileNotFoundError):
        sp = SPARCDocParser(".LaTeX")

    # Case 3: missing manual
    fs.create_file(".LaTeX/Introduction.tex")
    with pytest.raises(Exception):
        sp = SPARCDocParser(".LaTeX")

def test_docparser_init_working():
    """Mimic a working doc parser
    """
    from sparc.docparser import SPARCDocParser
    # Should work
    assert test_doc_dir.is_dir()
    sp = SPARCDocParser(directory=test_doc_dir)

    assert all([f.name != "SQ.tex" for f in sp.include_files])
    # No source code, date version will be None
    assert sp.version is None
    assert len(sp.parameters) > 0
    assert hasattr(sp, "parameter_categories")
    assert hasattr(sp, "other_parameters")
    # data type is a collected property only in to_json
    assert not hasattr(sp, "data_types")

def test_version_parser(fs, monkeypatch):
    """Only parse version"""
    from sparc.docparser import SPARCDocParser
    content_init_c = """void write_output_init(SPARC_OBJ *pSPARC) {
    int i, j, nproc, count;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // time_t current_time = time(NULL);
    // char *c_time_str = ctime(&current_time);
    time_t current_time;
    time(&current_time);
    char *c_time_str = ctime(&current_time);
    // ctime includes a newline char '\n', remove manually
    if (c_time_str[strlen(c_time_str)-1] == '\n') 
        c_time_str[strlen(c_time_str)-1] = '\0';

    FILE *output_fp = fopen(pSPARC->OutFilename,"w");
    if (output_fp == NULL) {
        printf("\nCannot open file \"%s\"\n",pSPARC->OutFilename);
        exit(EXIT_FAILURE);
    }

    fprintf(output_fp,"***************************************************************************\n");
    fprintf(output_fp,"*                       SPARC (version Feb 23, 2022)                      *\n");
    fprintf(output_fp,"*   Copyright (c) 2020 Material Physics & Mechanics Group, Georgia Tech   *\n");
    fprintf(output_fp,"*           Distributed under GNU General Public License 3 (GPL)          *\n");
    fprintf(output_fp,"*                   Start time: %s                  *\n",c_time_str);
    fprintf(output_fp,"***************************************************************************\n");
    fprintf(output_fp,"                           Input parameters                                \n");
    fprintf(output_fp,"***************************************************************************\n");
    if (pSPARC->Flag_latvec_scale == 0) {
        fprintf(output_fp,"CELL: %.15g %.15g %.15g \n",pSPARC->range_x,pSPARC->range_y,pSPARC->range_z);
        if (pSPARC->cell_typ <= 20) {
            fprintf(output_fp,"LATVEC:\n");
            fprintf(output_fp,"%.15f %.15f %.15f \n",pSPARC->LatUVec[0],pSPARC->LatUVec[1],pSPARC->LatUVec[2]);
            fprintf(output_fp,"%.15f %.15f %.15f \n",pSPARC->LatUVec[3],pSPARC->LatUVec[4],pSPARC->LatUVec[5]);
            fprintf(output_fp,"%.15f %.15f %.15f \n",pSPARC->LatUVec[6],pSPARC->LatUVec[7],pSPARC->LatUVec[8]);
    """
    
    fs.create_dir("src")
    fs.create_file("src/initialization.c")
    fs.create_dir("doc/.LaTeX")
    
    def mock_init(self):
        self.root = Path("doc/.LaTeX")
    monkeypatch.setattr(SPARCDocParser, "__init__", mock_init)
    sp = SPARCDocParser()
    with open("src/initialization.c", "w") as fd:
        fd.write(content_init_c)
    sp.parse_version(parse=False)
    assert sp.version is None
    sp.parse_version(parse=True)
    assert sp.version == "2022.02.23"

def test_include_files():
    """Test only include files"""
    from sparc.docparser import SPARCDocParser
    sp = SPARCDocParser(test_doc_dir)
    with pytest.warns(UserWarning):
        sp.get_include_files()
