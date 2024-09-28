# Folder for JOSS submission

This folder contains the `paper.md`, `paper.bib` and supporting
figures for the JOSS submission of the SPARC-X-API. Please do not add
SPARC doc in this directory.

## Using Building Actions



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
