# dmsa-pred

A standalone package for predicting variant phenotypes given a multiple sequence alignment and DMS data.

If this is used as a submodule to a nextflow build, use
we reccomend following [this guide](https://gist.github.com/gitaarik/8735255) as a workflow

When intergrating into an existing Nextstrain build, the primary steps to take are:

1. Add this as a submodule
2. Add the appropriate dependency of the resulting JSON's into the main workflow
3. install the dependencies:
```
$ conda activate ~/.nextstrain/runtimes/conda/env
$ snakemake --use-conda --cores 2 --create-envs-only --configfile dmsa-pipeline/dmsa-pred-reference.yaml
```
Now you can run the pipeline:
```
snakemake --use-conda --cores 4 --configfile dmsa-pipeline/dmsa-pred-reference.yaml
```




