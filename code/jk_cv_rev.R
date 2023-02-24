library(here)
library(revolver)
library(docopt)

'REVOLVER Analysis Script - Min Zhang, adapted by Jan Rogel

Usage:
  jk_cv_rev [options] 
  jk_cv_rev --version

Options:
  --min_group_size=<integer>  Minimum group size for the clusteirng algorithm to find. 
  --resamples=<integer>       Number of Jackknife resamplings.
  --version=1.0     Show version.

' -> doc

arguments <- docopt(doc)
# print(arguments)

results.dir    <- here("results")
jackknife.dir  <- here(results.dir,"Jackknife")
dir.create(jackknife.dir, recursive = T)

options.fit <- list(initial.solution = NA, 
                    transitive.orderings = FALSE, 
                    restarts = 10)

options.clustering <- list(transitive.closure = FALSE, 
                           min.group.size = 3, 
                           hc.method = 'ward', 
                           cutoff.features_annotation = 1, 
                           split.method = 'cutreeHybrid')

min.group.size <- as.integer(arguments$min_group_size)
resamples      <- as.integer(arguments$resamples)
print(min.group.size)
print(resamples)

# meso.fit <- readRDS("./input/meso_fitted_clone_trees.RDS")
meso.fit <- readRDS("../data/meso_fitted_clone_trees.RDS")

meso.fit.clustering <- revolver_cluster(
  meso.fit,
  split.method = 'cutreeHybrid',
  hc.method = "ward",
  min.group.size = min.group.size) # 3

# meso.fit.clustering$fit$clones_expansions$MED114
# meso.fit.clustering$fit$clones_expansions$MED150

# plot_trajectories_per_cluster(meso.fit.clustering, min_counts = 3)
# revolver::plot_clusters(meso.fit.clustering)


# plot_trajectories_per_cluster(meso.jackknife, min_counts = 3)

meso.jackknife <- revolver_jackknife(meso.fit.clustering, 
                                     resamples = resamples, # 100
                                     leave.out = 0.1, 
                                     options.fit = options.fit, 
                                     options.clustering = options.clustering)

saveRDS(meso.jackknife, here(jackknife.dir, paste0("meso.jackknife.",resamples,"rs.RDS")))


