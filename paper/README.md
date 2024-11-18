# Folder for JOSS submission

This folder contains the `paper.md`, `paper.bib` and supporting
figures for the JOSS submission of the SPARC-X-API. Please do not add
SPARC doc in this directory. While most of the writings should be done
in Markdown, please use either of the following methods for rendering
the draft into pdf format.

## Building Locally

Follow the instructions in the [JOSS
manual](https://joss.readthedocs.io/en/latest/paper.html#docker), use
docker (or equivalent) to build the paper locally:

```bash
# At the paper/ subfolder of the SPARC-X-API repository
docker run --rm \
    --volume $PWD:/data \
    --user $(id -u):$(id -g) \
    --env JOURNAL=joss \
    openjournals/inara
```
This command will create the `paper.pdf` under the `paper/` subfolder.

## Use Github Actions

The draft pdf will be rendered after any changes are pushed the
`paper/` under the `joss_paper` branch. Please check the [status
page](https://github.com/SPARC-X/SPARC-X-API/actions/workflows/joss_paper.yml)
for the latest build action. Once the compilation is finished, the zip
archive containing the rendered `paper.pdf` can be downloaded via the
link in the "Artifacts" section of the action status.
