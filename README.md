# Mesothelioma_evolution_deciphering_drugable_somatic_alterations

Analysis scripts and relative data used in the paper "Clonal Architecture in Mesothelioma is prognostic and shapes the Tumour Microenvironment" by Zhang et., al., 2022. https://www.nature.com/articles/s41467-021-21798-w

## data
Contains the input file (txt), which is already formatted and ready to use for REVOLVER analysis.

## R
This script performs REVOLVER analysis to construct phylogeny trees and infer repeat evolution trajectories in MEDUSA22 cohort.

```
Usage:
  revolver_analysis.R [options]

Options:
  --dataset=<FILE_PATH>     This is a required argument. Specify the path to the input file in .txt format - REQUIRED ARG
  --results=<FILE_PATH>     Speficy the outfolder to save the results in. - REQUIRED ARG
  --sspace=<VALUE>          Set the value for sspace.cutoff argument [default: 1000]
  --n_sampling=<VALUE>      [default: 500]
  --store_max=<VALUE>       [default: 200]
  --restarts=<VALUE>        [default: 10]
  --min_group_size=<VALUE>  [default: 3]
  --hc_method=<VALUE>       [default: ward]
  --split_method=<VALUE>    [default: cutreeHybrid]
  --tree_type               [defult:  clone_trees]
  -c --cohort               If set, looks for a constructed cohort in the ./ or results directory.
  -t --trees                If set, looks for computed trees in the ./ or results directory.
  -f --fit                  If set, looks for computed AND fitted trees in the ./ or results directory.
  -l --cluster              If set, looks for computed cohort, computed and fitted trees AND computed clusters. 
  -j --jackknife            If set, runs a quick non-intensive jackknife cross-validation.
  --version=0.1     Show version.

  For a basic run use (add -j to run cross-validation):

  Rscript revolver_analysis.R --dataset=input.txt --results=results -j 
```
