## Summary of page

The goal of this page is to identify variants of SARS-CoV-2 spike that are predicted to escape individual monoclonal antibodies or other therapeutics of interest.
We base our predictions on deep mutational scanning (DMS) experiments that measure the effects of mutations on antibody escape in the background of a single spike sequence (e.g., BA.1), using these effects to extrapolate predicted escape scores for all other spike sequences in the tree.
Each antibody of interest has an entry in "Color By" drop-down menu.
Selecting an entry will color the tree by predicted escape from that antibody.
The below sections list the antibodies we have analyzed, the sources of DMS data, and the strategy we use to predict escape.

## Approach for sampling SARS-CoV-2 sequences shown in the tree

The tree only shows a subset of available sequences, as the visualization can only handle ~4,000 sequences.
To downsample to a tractable number, we use the same sampling method as the [global ncov page](https://nextstrain.org/ncov/open/global/6m), which involves sampling sequences across the globe, enriching for sequences from the past 6 months.
Currently, we periodically update this page manually.
Though this will soon be automated to ensure that it is up to date with surveillance data.

## Antibodies of interest

We have organized antibodies into different datasets, accessible via the fourth drop-down menu under "Dataset". For a given dataset, the "Color By" drop-down menu has a separate entry for each antibody in the dataset.

The dataset called `unescaped-antibodies` includes all antibodies that we have analyzed that are likely to have appreciable activity against currently circulating variants.
Here is a list of these antibodies/therapeutics, along with the publication with DMS data used to predict escape:
* CC67.105 ([Dadonaite, 2023, Cell](https://doi.org/10.1016/j.cell.2023.02.001))
* CC9.104 ([Dadonaite, 2023, Cell](https://doi.org/10.1016/j.cell.2023.02.001))
* NTD_5-7 ([Dadonaite, 2023, Cell](https://doi.org/10.1016/j.cell.2023.02.001))
* C68.3 ([Guenthoer, 2022, biorxiv](https://doi.org/10.1101/2022.12.15.520606))
* C68.59 ([Guenthoer, 2022, biorxiv](https://doi.org/10.1101/2022.12.15.520606))
* C68.61 ([Guenthoer, 2022, biorxiv](https://doi.org/10.1101/2022.12.15.520606))
* LCB1 ([Hunt, 2022, Science Translational Medicine](https://doi.org/10.1126/scitranslmed.abn1252))

The other datasets group antibodies by publication, and include antibodies that the virus has already escaped.
Here is a list of these datasets:
* [Dadonaite, 2023, Cell](https://doi.org/10.1016/j.cell.2023.02.001)
* [Guenthoer, 2022, biorxiv](https://doi.org/10.1101/2022.12.15.520606)
* [Hunt, 2022, Science Translational Medicine](https://doi.org/10.1126/scitranslmed.abn1252)
* [Starr, 2022, PLoS Pathogens](https://doi.org/10.1371/journal.ppat.1010951)

## Naming of entries in the "Color By" drop-down menu

The options in the "Color By" drop-down menu are named by the antibody (prefix) followed by the type of predicted escape score (suffix).
There are two types of predicted escape scores, described in more detail below.
Their corresponding suffixes are:
* `_IC90_log_fold_change`: the predicted log<sub>10</sub> fold change in a sequence's IC<sub>90</sub> for a given antibody relative to the IC<sub>90</sub> of the DMS wildtype sequence.
* `_escape_score`: an escape score obtained by summing the effects of mutations that separate a given sequence from the DMS wildtype sequence. The effects are derived from the yeast-display DMS experiment described below.

For instance, `CC67.105_IC90_log_fold_change` corresponds to predictions for the C67.105 antibody, using the method for predicting changes in IC<sub>90</sub>.

Note: `C68.3-BA1` and `C68.3-BA2` are the same antibody with DMS escape measured in the BA.1 and BA.2 backgrounds, respectively.

## Strategies for predicting escape scores

Our strategy for making predictions is as follows.
For a given antibody, we first obtain a DMS dataset measuring the effects of mutations on antibody escape in the background of a given spike sequence.
Next, for each sequence in the tree (including inferred sequences of internal nodes), we make a list of amino-acid mutations that separate that sequence from the DMS sequence.
We then aggregate the effects of each mutation in the list, producing a single predicted escape score for that sequence.
If a mutation in the list does not have DMS data, then we ignore that mutation.
If a mutation in the list is at a site that does not align to a homologous site in Wuhan-Hu-1 spike (e.g., due to an insertion), then that mutation is also ignored (this is a limitation of our pipeline that we hope to address in the future).

The method we use to aggregate mutational effects depends on the DMS experiment.
We use one of two methods, listed below by naming suffix (see above):

1. `_IC90_log_fold_change`: The first method uses `polyclonal` ([Yu et al., 2022, Virus Evolution](https://doi.org/10.1093/ve/veac110)), which is a biophysical model that can be used to estimate mutational effects from a DMS experiment.
Given a set of inferred effects, `polyclonal` can also be used to estimate the IC<sub>90</sub> of unseen sequences for a given antibody.
We use this feature to estimate the IC<sub>90</sub> of each sequence in the tree.
We then compute the log<sub>10</sub> fold change in a sequence's IC<sub>90</sub> relative to the IC<sub>90</sub> of the DMS sequence, and use this metric to color trees.
The corresponding color scale always ranges from -1 to 4 to provide a consistent scale across antibodies, though values reported in the tool tip can exceed this range.

2. `_escape_score`: The second method involves adding mutational effects quantified using a yeast-display DMS experiment devised by [Greaney et al., 2021, Cell Host & Microbe](https://doi.org/10.1016/j.chom.2020.11.007).
In this experiment, mutant variants of the RBD are displayed on the surface of yeast, with each cell expressing many copies of a single variant.
Cells are then labeled with an anti-RBD antibody and fluorescence-activated cell sorting is used to select cells with RBD variants that escape antibody binding.
Deep sequencing is then used to compute the frequency of each mutation in the pool of cells before selection and the antibody-escaped pool of cells after selection.
The resulting data are used to compute an "escape fraction" for each mutation.
We aggregate these mutational effects by simply summing the escape fractions of individual mutations.
When coloring trees with this metric, the color scale always ranges from 0 to 1 to provide a consistent scale across antibodies, though values reported in the tool tip can exceed this range.
