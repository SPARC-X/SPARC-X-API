# Setup Environment for SPARC-X-API

`SPARC-X-API` is designed to automate the discovery of
pseudopotential files, the JSON API, and the SPARC binary. However,
you can exert fine-grained control over their setup:

### A) Pseudopotential files
Pseudopotential files (in `Abinit` psp8 format) are loaded in the following
order:

1) Via the `psp_dir` argument passed to the `sparc.SPARC` calculator.
2) Through the environment variables `$SPARC_PSP_PATH` or `$SPARC_PP_PATH` (this is the
 method employed by [`conda` installation](#1-via-anaconda-or-miniconda-recommended)).
3) By using `psp8` files bundled with the SPARC-X-API installation (see the
[manual installation](#2-manual-installation-from-source-with-pip)).

To specify a custom path for your psp8 files, set the `$SPARC_PSP_PATH` or `$SPARC_PP_PATH` variable as follows:
```bash
export SPARC_PSP_PATH="/path/to/your/psp8/directory"
```

To determine the default location of psp8 files (as per option 3), run the following code:
```bash
python -c "from sparc.common import psp_dir; print(psp_dir)"
```

### B) JSON schema
`SPARC-X-API` is engineered for compatibility with the SPARC
C-code. It achieves this by loading a JSON schema for
parameter validation and unit conversion. You can review the default
schema used by the API at sparc.sparc_json_api.default_json_api

```json
"FD_GRID": {
   "symbol": "FD_GRID",
   "label": "FD_GRID",
   "type": "integer array",
   "default": null,
   "unit": "No unit",
   "example": "FD_GRID: 26 26 30",
   "description": "#<Some description...>",
   "allow_bool_input": false,
   "category": "system"
  },
```

The schema file is generated from SPARC's LaTeX documentation.  In
upcoming releases of `SPARC-X-API`, we're aiming to provide users the
flexibility to use their own custom schema files. This would be
particularly useful for those who might be testing a development
branch of SPARC. By default, the JSON schema is packaged under
`sparc/sparc_json_api` directory. If you have another version of SPARC
source code, you can set the environment variable `$SPARC_DOC_PATH` to
the directory containing the LaTeX codes for the documentation, such
as `<SPARC-source-code-root>/doc/.LaTeX`. If you obtain `sparc-x` from
the conda method as mentioned above, By default, the JSON schema is
packaged under `sparc/sparc_json_api` directory. If you have another
version of SPARC source code, you can set the environment variable
`$SPARC_DOC_PATH` is automatically set to
`<conda-env-root>/share/doc/sparc/.LaTeX`. Setting up the environment
variable `$SPARC_DOC_PATH` helps loading the correct JSON schame that
is compatible with your SPARC binary code.

### C) SPARC Command Configuration

The command to execute SPARC calculations is determined based on the following priority:

1) The command argument provided directly to the `sparc.SPARC` calculator.
2) The environment variable `$ASE_SPARC_COMMAND`
3) If neither of the above is defined, `SPARC-X-API` looks for the SPARC binary under current `$PATH` and combine with the suitable `mpi` command prefix.

Example:

1. Using `mpirun` (e.g. on a single test machine)
```bash
export ASE_SPARC_COMMAND="mpirun -n 8 -mca orte_abort_on_non_zero_status 1 /path/to/sparc -name PREFIX"
```

2. Using `srun` (e.g. in HPC slurm job system)
```bash
export ASE_SPARC_COMMAND="srun -n 8 --kill-on-bad-exit /path/to/sparc -name PREFIX"
```

*Notes*:
1. The `-name PREFIX` part can be omitted the `label` property of the `sparc.SPARC` calculator is set (which is the default behavior). Any extra features of the SPARC code (e.g. GPU acceleration) should be specified in the command.

2. We recommend adding kill switches for your MPI commands like the examples above when running `sparc` to avoid unexpected behaviors with exit signals.
