# dmsa-pred

A standalone package for predicting variant phenotypes given a multiple sequence alignment and DMS data.

If this is used as a submodule to a nextflow build, use
```
git submodule update --init --recursive
```
To update the submodule.

When intergrating into an existing Nextstrain build, the primary steps to take are:

1. Add this as a submodule
2. Add the appropriate dependency of the resulting JSON's into the main workflow
3. install the dependencies:
```
$ nextstrain shell
$ snakemake --use-conda --cores 2 --create-envs-only --configfile dmsa-pipeline/dmsa-pred-reference.yaml
```

