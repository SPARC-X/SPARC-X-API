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

You can configure the setup for SPARC-X-API using either of the
following methods:
1. (Recommended) use the [ASE configuration file](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html#calculator-configuration)
2. Use environmental variables. Please note although environmental
   variables have long been the standard approach to set
   ASE-calculators, they may be obsolete in future releases of ASE.

```{note}
The environmental variables will have **higher priority** than the equivalent fields in the configure file, if both are set.
```

### Editing the configuration file

ASE will look for a configuration file at `~/.config/ase/config.ini`
for package-specific settings. The configuration file follows the [INI
format](https://en.wikipedia.org/wiki/INI_file), where key-value pairs
are grouped in sections. An example of the SPARC-specific section may look like follows:

```{code} ini
[sparc]
; `command`: full shell command (include MPI directives) to run SPARC calculation
;            has the same effect as `ASE_SPARC_COMMAND`
command = srun -n 24 ~/bin/sparc

; `psp_path`: directory containing pseudopotential files
;            has the same effect as `SPARC_PSP_PATH`
psp_path = ~/dev_SPARC/psps

; `doc_path`: directory for SPARC LaTeX documentation
;             has the same effect as `SPARC_DOC_PATH`
doc_path = ~/dev_SPARC/doc/.LaTeX/
```

The available options in the configuration file are:
1. SPARC command: either use `command` to set a full shell string to
   run the SPARC program, or use the combination of `sparc_exe` and
   `mpi_prefix`. See [SPARC command configuration](#sparc-cmd-setting)
   for more details.
2. JSON schema settings (*Optional*): use either `json_schema` to
   define a custom JSON schema file, or `doc_path` for parsing the
   LaTeX documentation on-the-fly. See [JSON schema
   configuration](#json-schema-setting) for more details.
3. Pseudopotential settings (*Optional*): use `psp_path` for the
   location of pseudopotential files. See [pseudopotential files settings](#pseudopot-setting)
   for more details.


You can overwrite the location of the configuration file by the
environmental variable `ASE_CONFIG_PATH`.

(json-schema-setting)=
### JSON schema

Each version of SPARC-X-API ships with a JSON schema compatible with a
dated-version of SPARC C/C++ code. You can find that version in the README badge like follows:

![API badge](https://raw.githubusercontent.com/SPARC-X/SPARC-X-API/badges/badges/api_version.svg)

If that does not match your local SPARC version, you can configure the location of JSON schema using one of the following methods:

#### Option 1. Parse LaTeX documentation on-the-fly

The environment variable `SPARC_DOC_PATH` (equivalent to `doc_path` field in configuration file)  will direct SPARC-X-API to look for a local directory containing LaTeX documentation to parse on-the-fly, for example:
```bash
export SPARC_DOC_PATH=<local-SPARC-dir>/doc/.LaTeX
```
or configuration file setting:
```{code} ini
doc_path: <local-SPARC-dir>/doc/.LaTeX
```


#### 2. Use your own `parameters.json`

In some cases an experimental feature may not have been updated in the
official doc. You can create and edit your own `parameters.json` file
to temporarily test a local version of SPARC:

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

Finally, set the `json_schema` field in the configuration file to the
newly generated json file, for example:
```{code} ini
[sparc]
; `json_schema`: custom schema file parsed from LaTeX documentation
json_schema: ~/SPARC/parameters.json
```

```{warning}
`json_schema` and `doc_path` fields cannot be both present in the configuration file!
```

(pseudopot-setting)=
### Pseudopotential files


To specify a custom path for your pseudopotential files (in Abinit [psp8 format](https://docs.abinit.org/developers/psp8_info/#:~:text=The%20format%208%20for%20ABINIT,great%20flexibility%20in%20doing%20so.),
you can either use the environment variable `SPARC_PSP_PATH`:
```bash
export SPARC_PSP_PATH=path/to/your/psp8/directory
```
or the equivalent keyword `psp_path` in configuration file
```{code} ini
psp_path: path/to/your/psp8/directory
```

When installing SPARC [via `conda-forge`](#use-conda),
`$SPARC_PSP_PATH` is already included in the activate script of the
conda environment.


To determine the location of default psp8 files (as in [manual pip installation](#pip-install)), run the following code:
```bash
python -c "from sparc.common import psp_dir; print(psp_dir)"
```

(sparc-cmd-setting)=
### SPARC command configuration

The command to execute SPARC calculations is determined based on the
following priority:

1. The command argument provided directly to the `sparc.SPARC` calculator.
2. Command defined by environment variable `ASE_SPARC_COMMAND` or configuration file
3. If neither of the above is defined, `SPARC-X-API` looks for the SPARC binary under current `$PATH` and combine with the suitable MPI command prefix.

#### Use full command

The variable `ASE_SPARC_COMMAND` (or `command` field in configuration file)
contain the *full command* to run a SPARC calculation. depending on the system, you may choose one of the following:

1. Using `mpirun` (e.g. on a single test machine)
```bash
export ASE_SPARC_COMMAND="mpirun -n 8 -mca orte_abort_on_non_zero_status 1 /path/to/sparc -name PREFIX"
```
or in configuration file
```{code} ini
command: mpirun -n 8 -mca orte_abort_on_non_zero_status 1 /path/to/sparc -name PREFIX
```

2. Using `srun` (e.g. [SLURM](https://slurm.schedmd.com/documentation.html) job system HPCs)
```bash
export ASE_SPARC_COMMAND="srun -n 8 --kill-on-bad-exit /path/to/sparc -name PREFIX"
```
or in configuration file
```{code} ini
command: srun -n 8 --kill-on-bad-exit /path/to/sparc -name PREFIX" /path/to/sparc -name PREFIX
```

```{note}
1. The `-name PREFIX` is optional and will automatically replaced by the `sparc.SPARC` calculator.

2. We recommend adding kill switches for your MPI commands like the examples above when running `sparc` to avoid unexpected behaviors with exit signals.
```

#### Specifying MPI binary location

It is also possible to construct the SPARC command from two

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
2. Python API documentation: https://github.com/SPARC-X/SPARC-X-API/blob/master/README.md

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
