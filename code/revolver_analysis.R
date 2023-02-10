###########################################################
##                                                       ##
## Package: revolver                                     ##
## Version: 0.2.0                                        ##
## Date: 2013-10-2                                       ##
## Title: Repeated Evolution in Cancer (REVOLVER-updates)##
##                                                       ##
###########################################################
library(here)
library(revolver)

## options 
options(crayon.enabled = FALSE)
options(revolver.progressBar = FALSE)
options(revolver.jamPDF=TRUE)

options.trees = list(sspace.cutoff = 1000, n.sampling = 500, store.max = 200, overwrite = FALSE)
options.fit = list(initial.solution = NA, transitive.orderings = FALSE, restarts = 10)
options.clustering.withGL = list(transitive.closure = FALSE, min.group.size = 3, hc.method = 'ward', cutoff.features_annotation = 1, split.method = 'cutreeHybrid')
options.clustering.withoutGL = list(transitive.closure = FALSE, min.group.size = 3, hc.method = 'ward', cutoff.features_annotation = 1,split.method = 'cutreeDynamic')

load(here("data", "revolver_analysis.Rdata"))

cohort.name = 'MEDUSA'
do.plots = TRUE
cex.fit = 3
folder.output = "./"
phylo.dir = file.path(folder.output,"Phylogenies")
clustering.dir = file.path(folder.output,"Clustering")
jackknife.dir = file.path(folder.output,"Jackknife")

dir.create(phylo.dir, showWarnings = F, recursive = T)
dir.create(clustering.dir, showWarnings = F, recursive = T)
dir.create(jackknife.dir, showWarnings = F, recursive = T)

## construct revolver object
meso.cohort = revolver_cohort(
  dataset = dataset,
  CCF_parser = revolver:::CCF_parser,
  annotation = cohort.name,
  ONLY.DRIVER = FALSE,
  MIN.CLUSTER.SIZE = 1
)

#----- check cohort for singular drivers and remove them ----
removeSingularDrivers <- function(meso.cohort, stopOnError = TRUE){
  non.recurrent = Stats_drivers(meso.cohort) %>%
    filter(N_tot == 1) %>%
    pull(variantID)
  
  for(VID in non.recurrent){
    meso.cohort <- remove_drivers(meso.cohort, variantID = VID, check = TRUE)
  }
  
  revolver_check_cohort(meso.cohort, stopOnError = TRUE)
  
  return(meso.cohort)
}
meso.cohort <- removeSingularDrivers(meso.cohort = meso.cohort)
revolver_check_cohort(meso.cohort, stopOnError = TRUE)
non.recurrent = Stats_drivers(meso.cohort) %>%
  filter(N_tot == 1) %>%
  pull(variantID)

for(VID in non.recurrent){
  meso.cohort <- remove_drivers(meso.cohort, variantID = VID, check = TRUE)
}
revolver_check_cohort(meso.cohort, stopOnError = TRUE)


## infer phylogeny tree of subclonal/clonal
for (patient in meso.cohort$patients) {
  meso.cohort <- compute_mutation_trees(meso.cohort,patient, options.trees)
}

meso.cohort <- compute_mutation_trees(meso.cohort, 
                                      patients     = meso.cohort$patients, 
                                      sspace.cutof = 1000,
                                      n.sampling   = 500,
                                      store.max    = 200, 
                                      overwrite    = FALSE)

## fit repeat evolution tree
meso.fit <- revolver_fit(meso.cohort, 
                        initial.solution = options.fit$initial.solution)
# --- cluster the tree
meso.fit.clustering <- revolver_cluster(
  meso.fit, 
  split.method = 'cutreeHybrid',
  hc.method = "ward",
  min.group.size = 3)

plot_DET_index(meso.fit)
plot_clusters(meso.fit.clustering)
plot_drivers_clonality(meso.fit)
plot_drivers_graph(meso.fit)
plot_drivers_occurrence(meso.fit)
plot_penalty(meso.fit)
plot_jackknife_cluster_stability(meso.jackknife)
plot_jackknife_coclustering(meso.jackknife)
plot_jackknife_trajectories_stability(meso.jackknife)



# for (patient in meso.cohort$patients) {
# 	patient.dir = file.path(phylo.dir,patient)
# 	dir.create(patient.dir, showWarnings = F, recursive = T)
# 	f = file.path(patient.dir,paste(patient, c("fit.pdf", "trajectories.pdf", "itransfer.pdf"),sep="-"))
# 	revolver_plt_fit_patient(meso.fit, patient, cex = cex.fit, file = f[1])
# 	revolver_plt_trajectories_patient(meso.fit, patient, cex = cex.fit, file = f[2])
# 	revolver_plt_itransfer_patient(meso.fit, patient, cex = cex.fit, file = f[3])
# }

# meso.fit.clustering = revolver_evo_distance(meso.fit, use.GL = TRUE, transitive.closure = options.clustering.withGL$transitive.closure)
# revolver_infoclustering(meso.fit.clustering, min.group.size = options.clustering.withGL$min.group.size, 
#     do.plot = do.plots, file = file.path(clustering.dir,"Cluster_gl-infoClusterting-report.pdf"))
# meso.fit.clustering = revolver_cluster(meso.fit.clustering, hc.method = options.clustering.withGL$hc.method, 
#     min.group.size = options.clustering.withGL$min.group.size, 
#     split.method = options.clustering.withGL$split.method)

## plot cluster and heatmap
# pdf(file.path(clustering.dir,"plt_rclusters.GL.pdf"),width=18,height=8)
# revolver_plt_rclusters(meso.fit.clustering,cutoff.features_annotation = 1,symbol.arrow = "->")
# dev.off()
# pdf(file.path(clustering.dir,"plt_dendrogram.GL.pdf"),width=6,height=4)
# revolver_plt_dendrogram(meso.fit.clustering)
# dev.off()
# revolver_plt_evodistance(meso.fit.clustering,file=file.path(clustering.dir,"plt_evodistance.GL.pdf"),cex=3)

meso.fit.clustering

## jackknife 
meso.jackknife = revolver_jackknife(meso.fit.clustering)
# meso.jackknife = revolver_jackknife(meso.fit.clustering) 
Stats_fits()
Stats_drivers(meso.cohort)
write.csv(Stats_drivers(meso.cohort),"cohort's_driver_events.csv")
Stats_fits(meso.fit)
write.csv(Stats_fits(meso.cohort),"cohort's_fits.csv")
Stats_trees(meso.cohort)
write.csv(Stats_trees(meso.cohort),"cohort's trees.csv")
DET_index(meso.fit)
write.csv(DET_index(meso.cohort),"index_of_DET.csv")
Jackknife_patient_coclustering(meso.jackknife)
write.csv(Jackknife_patient_coclustering(meso.cohort),"jackknife_patient_coclustering.csv")
Jackknife_cluster_stability(meso.jackknife)
write.csv(Jackknife_cluster_stability(meso.cohort),"jackknife_cluster_stability.csv")
Jackknife_trajectories_stability(meso.jackknife)
write.csv(Jackknife_trajectories_stability(meso.cohort),"jackknife_trajectories_stability.csv")

revolver_plt_jackknife_coclust(
	x = meso.jackknife,
	cutoff.annotate.numbers = 0.6,
	file = file.path(jackknife.dir,"MEDUSA.jackknife_coclust.pdf")
)
revolver::plot
revolver_plt_jackknife_coclust_bplot(
	x = meso.jackknife,
	file = file.path(jackknife.dir,"MEDUSA.jackknife_coclust_bplot.pdf")
)

revolver_plt_jackknife_edge_counts(
	x = meso.jackknife,
	cutoff.annotate.numEdges = 3,
	file = file.path(jackknife.dir,"MEDUSA.jackknife_edge_counts.pdf")
)

revolver_plt_jackknife_edge_prb(
	x = meso.jackknife,
	sorting = "heterogeneity",
	cutoff.annotate.numbers = 0.6,
	file = file.path(jackknife.dir,"MEDUSA.jackknife_edge_prb.pdf")
)
