# Configurations for SPARC-X-API

`SPARC-X-API` requires the following components to be fully functional:
- A JSON schema parsed from SPARC LaTeX documentation
- Pseudopotential files
- SPARC C/C++ binary

## Default configurations

`SPARC-X-API` is designed to automate the discovery of these
configurations. The default configurations are:

- JSON schema: `<sparc-x-api-root>/sparc_json_api/parameters.json`

- Pseudopotential files: `<sparc-x-api-root>/psp/*.psp8` (if installed [via conda](#use-conda))

- Command to run SPARC binary: `mpirun sparc` (if `sparc` is found in `$PATH`)

## Custom configurations

**TODO** use sparc ini config.

Please check the following steps for detailed setup:

### JSON schema

Each version of SPARC-X-API ships with a JSON schema compatible with a
dated-version of SPARC C/C++ code. You can find that version in the README badge like follows:

![API badge](https://raw.githubusercontent.com/SPARC-X/SPARC-X-API/badges/badges/api_version.svg)

If that does not match your local SPARC version, you can configure the location of JSON schema using one of the following methods:

1. Set the `$SPARC_DOC_PATH` variable

`$SPARC_DOC_PATH` will direct SPARC-X-API to look for a local directory containing LaTeX documentation, for example:
```bash
export SPARC_DOC_PATH=<local-SPARC-dir>/doc/.LaTeX
```
and parse the JSON schema on-the-fly.

2. Use your own `parameters.json`

In some cases an experimental feature may not have been updated in the
official doc. You can create and edit your own `parameters.json` file to
temporarily test a local version of SPARC:

First parse the LaTeX files into `parameters.json`
```bash
python -m sparc.docparser  <local-SPARC-dir>/doc/.LaTeX \
                           --output parameters.json \
                           --include-subdirs
```

Then add / edit missing parameters in the `parameters` section in
`parameters.json`, see examples from the existing file:
```json
"ACE_FLAG": {
   "symbol": "ACE_FLAG",
   "label": "ACE_FLAG",
   "type": "integer",
   "default": 1,
   "unit": "No unit",
   "example": "ACE_FLAG: 0",
   "description": "Use ACE operator to accelarte the hybrid calculation.",
   "remark": "Without ACE operator, the hybrid calculation will be way slower than with it on depending on the system size.",
   "allow_bool_input": true,
   "default_remark": "1",
   "description_raw": "Use ACE operator to accelarte the hybrid calculation.",
   "remark_raw": "Without ACE operator, the hybrid calculation will be way slower than with it on depending on the system size.",
   "category": "scf"
  },
```

**TODO** make sure the json parameters are included.
Finally, set up **TODO**.

<!--t `SPARC-X-API` is engineered for compatibility with the SPARC -->
<!-- C-code. It achieves this by loading a JSON schema for --> <!--
parameter validation and unit conversion. You can review the default
--> <!-- schema used by the API at
sparc.sparc_json_api.default_json_api -->

<!-- ```json -->
<!-- "FD_GRID": { -\-> -->
<!-- "symbol": "FD_GRID", -\-> -->
<!-- <\!--    "label": "FD_GRID", -\-> -->
<!-- <\!--    "type": "integer array", -\-> -->
<!-- <\!--    "default": null, -\-> -->
<!-- <\!--    "unit": "No unit", -\-> -->
<!-- <\!--    "example": "FD_GRID: 26 26 30", -\-> -->
<!-- <\!--    "description": "#<Some description...>", -\-> -->
<!-- <\!--    "allow_bool_input": false, -\-> -->
<!-- <\!--    "category": "system" -\-> -->
<!-- <\!--   }, -\-> -->
<!-- <\!-- ``` -\-> -->

<!-- <\!-- The schema file is generated from SPARC's LaTeX documentation.  In -\-> -->
<!-- <\!-- upcoming releases of `SPARC-X-API`, we're aiming to provide users the -\-> -->
<!-- <\!-- flexibility to use their own custom schema files. This would be -\-> -->
<!-- <\!-- particularly useful for those who might be testing a development -\-> -->
<!-- <\!-- branch of SPARC. By default, the JSON schema is packaged under -\-> -->
<!-- <\!-- `sparc/sparc_json_api` directory. If you have another version of SPARC -\-> -->
<!-- source code, you can set the environment variable `$SPARC_DOC_PATH` to -->
<!-- the directory containing the LaTeX codes for the documentation, such -->
<!-- as `<SPARC-source-code-root>/doc/.LaTeX`. If you obtain `sparc-x` from -->
<!-- the conda method as mentioned above, By default, the JSON schema is -->
<!-- packaged under `sparc/sparc_json_api` directory. If you have another -->
<!-- version of SPARC source code, you can set the environment variable -->
<!-- `$SPARC_DOC_PATH` is automatically set to -->
<!-- `<conda-env-root>/share/doc/sparc/.LaTeX`. Setting up the environment -->
<!-- variable `$SPARC_DOC_PATH` helps loading the correct JSON schame that -->
<!-- is compatible with your SPARC binary code. -->


### Pseudopotential files

Pseudopotential files (in `Abinit` [psp8 format](https://docs.abinit.org/developers/psp8_info/) are loaded in the
following order:

<!-- 1. Via the `psp_dir` argument passed to the `sparc.SPARC` calculator. -->
<!-- 2. Through the environment variables `$SPARC_PSP_PATH` or `$SPARC_PP_PATH` (this is the -->
<!--  method employed by [`conda` installation](#use-conda)). -->
<!-- (manual-psp8)= -->
<!-- 3. By using `psp8` files bundled with the SPARC-X-API installation (see the -->
<!-- [pip installation](#pip-install)). -->

To specify a custom path for your pseudopotential files (in Abinit [psp8 format]()),
use the environment variable `$SPARC_PSP_PATH` or `$SPARC_PP_PATH` variable:
```bash
export SPARC_PSP_PATH=path/to/your/psp8/directory
```

**TODO** we should sunset the `SPARC_PP_PATH`

When installing SPARC [via `conda-forge`](#use-conda),
`$SPARC_PSP_PATH` is already included in the activate script of the
conda environment.


To determine the location of default psp8 files (as in [manual pip installation](#pip-install)), run the following code:
```bash
python -c "from sparc.common import psp_dir; print(psp_dir)"
```


### SPARC Command Configuration

The command to execute SPARC calculations is determined based on the
following priority:

1. The command argument provided directly to the `sparc.SPARC` calculator.
2. The environment variable `$ASE_SPARC_COMMAND`
3. If neither of the above is defined, `SPARC-X-API` looks for the SPARC binary under current `$PATH` and combine with the suitable MPI command prefix.

Example to set `$ASE_SPARC_COMMAND`

1. Using `mpirun` (e.g. on a single test machine)
```bash
export ASE_SPARC_COMMAND="mpirun -n 8 -mca orte_abort_on_non_zero_status 1 /path/to/sparc -name PREFIX"
```

2. Using `srun` (e.g. [SLURM](https://slurm.schedmd.com/documentation.html) job system HPCs)
```bash
export ASE_SPARC_COMMAND="srun -n 8 --kill-on-bad-exit /path/to/sparc -name PREFIX"
```

```{note}
1. The `-name PREFIX` is optional and will automatically replaced by the `sparc.SPARC` calculator.

2. We recommend adding kill switches for your MPI commands like the examples above when running `sparc` to avoid unexpected behaviors with exit signals.
```

## Post-installation check

We recommend the users to run a simple test after installation and
setup to make sure everything works:

```bash
python -m sparc.quicktest
```

A proper setup will display the following in the output's Summary section:
```{raw} html
<div class="highlight">
<pre>
--------------------------------------------------------------------------------
<font color="#AD7FA8"><b>Summary</b></font>
--------------------------------------------------------------------------------
<font color="#AD7FA8"><b>Configuration</b></font>
psp_dir: /home/pink/Dev/SPARC/psps
api_version: 2024.10.14
api_source: {&apos;path&apos;: &apos;/home/pink/Dev/SPARC-X-API/sparc/sparc_json_api/parameters.json&apos;, &apos;type&apos;: &apos;json&apos;}
command: mpirun -n 2 /home/pink/Dev/dev_SPARC/lib/sparc
sparc_version: 2024.10.14
sparc_socket_compatibility: True
--------------------------------------------------------------------------------
<font color="#AD7FA8"><b>Tests</b></font>
<b>Import:</b><font color="#8AE234">  PASS</font>
<b>Pseudopotential:</b><font color="#8AE234">  PASS</font>
<b>JSON API:</b><font color="#8AE234">  PASS</font>
<b>SPARC Command:</b><font color="#8AE234">  PASS</font>
<b>Calculation (File I/O):</b><font color="#8AE234">  PASS</font>
<b>Calculation (UNIX socket):</b><font color="#8AE234">  PASS</font>
--------------------------------------------------------------------------------
</pre>
</div>
```

Check for error messages when some tests didn't pass and
troubleshooting hints, such as the example below with a mis-configured command.

```{raw} html
<div class="highlight">
<pre>--------------------------------------------------------------------------------
<font color="#AD7FA8"><b>Summary</b></font>
--------------------------------------------------------------------------------
<font color="#AD7FA8"><b>Configuration</b></font>
psp_dir: /home/pink/Dev/SPARC/psps
api_version: 2024.10.14
api_source: {&apos;path&apos;: &apos;/home/pink/Dev/dev_SPARC/doc/.LaTeX&apos;, &apos;type&apos;: &apos;latex&apos;}
command: mpirun -n 4 /home/pink/Dev/dev_SPARC/lib/sparc
sparc_version: NaN
sparc_socket_compatibility: False
--------------------------------------------------------------------------------
<font color="#AD7FA8"><b>Tests</b></font>
<b>Import:</b><font color="#8AE234">  PASS</font>
<b>Pseudopotential:</b><font color="#8AE234">  PASS</font>
<b>JSON API:</b><font color="#8AE234">  PASS</font>
<b>SPARC Command:</b><font color="#EF2929">  FAIL</font>
<b>Calculation (File I/O):</b><font color="#EF2929">  FAIL</font>
<b>Calculation (UNIX socket):</b><font color="#EF2929">  FAIL</font>
--------------------------------------------------------------------------------
<font color="#EF2929">Some tests failed! Please check the following information.</font>

<b>SPARC Command:</b>
<font color="#EF2929">Error detecting SPARC version</font>
- The command prefix to run SPARC calculation should look like
  `&lt;mpi instructions&gt; &lt;sparc binary&gt;`
- Use $ASE_SPARC_COMMAND to set the command string
- Check HPC resources and compatibility (e.g. `srun` on a login node)


<b>Calculation (File I/O):</b>
<font color="#EF2929">Simple calculation in file I/O mode failed: </font>
<font color="#EF2929">SPARC failed with command mpirun -n 4 /home/pink/Dev/dev_SPARC/lib/sparc -name SPARCwith error code 1</font>
- Check if settings for pseudopotential files are correct
- Check if SPARC binary exists and functional
- Check if specific HPC requirements are met:
  (module files, libraries, parallel settings, resources)


<b>Calculation (UNIX socket):</b>
<font color="#EF2929">Simple calculation in socket mode (UNIX socket) failed: </font>
<font color="#EF2929">Cannot find the sparc executable! Please make sure you have the correct setup</font>
- The same as error handling in file I/O calculation test
- Check if SPARC binary supports socket


--------------------------------------------------------------------------------
Please check additional information from:
1. SPARC&apos;s documentation: https://github.com/SPARC-X/SPARC/blob/master/doc/Manual.pdf
2. Python API documentation: https://github.com/alchem0x2A/SPARC-X-API/blob/master/README.md

</pre>
</div>
```

```{note}
1. When using SPARC-X-API to parse SPARC files, it's essential that at
least the "Import" and "JSON API" tests are successful.
2. For running
SPARC calculations, "SPARC Command" and "Calculation (File I/O)" must also
succeed.
3. "Calculation (UNIX socket)" ensures the SPARC binary is compatible with socket communication,
see [calculation in socket mode](advanced_socket.md).
```

If you run into further problems, consult our [troubleshooting
guidlines](troubleshooting.md) or [raise an issue](contribute.md).
